#!/usr/bin/env nextflow

// Create queue channels (consumable)
// will need to split up into tumor channel and normal channel, use regex for this
tumor_ch = Channel.fromPath("${params.bam_files}/*_T_*.bam")
normal_ch = Channel.fromPath("${params.bam_files}/*_G_*.bam")


// Create value channels (use first operator to convert to value from queue)
tumor_val = tumor_ch.first()
normal_val = normal_ch.first()

// Sample id channel: take full filename and grab the sample id by splitting by _ and taking the first item
sample_id = tumor_ch.map { filePath -> 
                def fileName = filePath.baseName
                def baseName = fileName.split('_')[0]
                return baseName}
sample_id_ch = sample_id.first()

// Define the list of chromosomes + create a channel emitting each chromosome
chromosomes = (1..2).collect { it.toString() } + ['X']
chrom_strings = Channel.from(chromosomes)
chrom_ch = chrom_strings.map { it -> "chr" + it }


include { GETPILEUPSUMMARIES } from '../../tools/gatk/get_pileup_summaries.nf'
include { CALCULATECONTAMINATION } from '../../tools/gatk/calculate_contamination.nf'
include { MUTECT2 } from '../../tools/gatk/mutect.nf'
include { BGZIP; PREPAREVCF } from '../../tools/bcftools/prepareVCFs.nf'
include { MERGESTATS } from '../../tools/bcftools/combineMutectStats.nf'

workflow {
    /*
    GETPILEUPSUMMARIES(tumor_ch, normal_ch, params.exac)
    tumor_table = GETPILEUPSUMMARIES.out.tumor
    normal_table = GETPILEUPSUMMARIES.out.normal
    
    CALCULATECONTAMINATION(tumor_table, normal_table)
    contam_table = CALCULATECONTAMINATION.out.contamination
    segment_table = CALCULATECONTAMINATION.out.segment
     
    */
    // Run mutect2
    MUTECT2(tumor_val, normal_val, chrom_ch, sample_id_ch)
    
    // Merge and prepare VCF
    BGZIP(MUTECT2.out.vcf)
    vcfs_ch = BGZIP.out.vcf.collect()
    PREPAREVCF(vcfs_ch, sample_id_ch)
    
    // Merge stats, map the command line call in MergeMutectStats to each path
    stats = MUTECT2.out.stats.map { it -> "-stats " + it}
    stats_ch = stats.collect()
    MERGESTATS(stats_ch)
}