#!/usr/bin/env nextflow

// Create queue channels (consumable)
// will need to split up into tumor channel and normal channel, use regex for this

all_bams = Channel.fromPath("${params.bam_files}/*.bam")

// Define the list of chromosomes + create a channel emitting each chromosome
chromosomes = (1..22).collect { it.toString() } + ['X']
chrom_strings = Channel.from(chromosomes)
chrom_ch = chrom_strings.map { it -> "chr" + it }

include { RENAME } from '../../tools/samtools/get_samplename.nf'
include { SORT; SORTANDINDEX } from '../../tools/samtools/sort_and_index.nf'
include { MARKDUPLICATES } from '../../tools/gatk/mark_duplicates.nf'
include { GETPILEUPSUMMARIES } from '../../tools/gatk/get_pileup_summaries.nf'
include { CALCULATECONTAMINATION } from '../../tools/gatk/calculate_contamination.nf'
include { MUTECT2 } from '../../tools/gatk/mutect.nf'
include { BGZIP; PREPAREVCF } from '../../tools/bcftools/prepareVCFs.nf'
include { MERGESTATS } from '../../tools/bcftools/combineMutectStats.nf'
include { LEARNORIENTATION } from '../../tools/bcftools/combineF1R2files.nf'
include { FILTERMUTECT } from '../../tools/gatk/filter_mutect.nf'
include { PASS } from '../../tools/snpeff/sift_variants.nf'
include { ANNOTATE } from '../../tools/snpeff/annotate_variants.nf'

workflow {
    RENAME(all_bams)

    // sort with samtools 
    SORT(RENAME.out)

    // mark duplicates
    MARKDUPLICATES(SORT.out)

    // sort and index with samtools to prep for gatk somatic variant calling
    SORTANDINDEX(MARKDUPLICATES.out.bam)
    
    bam_dir = SORTANDINDEX.out.bam
    bai_dir = SORTANDINDEX.out.bai

    // gatk getpileupsummaries
    GETPILEUPSUMMARIES(bam_dir, params.exac)
    tumor_table = GETPILEUPSUMMARIES.out.tumor
    normal_table = GETPILEUPSUMMARIES.out.normal.first()
    
    // gatk calculate contamination from pileup summaries
    CALCULATECONTAMINATION(tumor_table, normal_table)
    contam_table = CALCULATECONTAMINATION.out.contamination.collect()
    segment_table = CALCULATECONTAMINATION.out.segment.collect()
    
    // separate out tumor and normal samples into two different channels
    def tumorpattern = params.tumor
    def normalpattern = params.normal
    tumor_ch = bam_dir.filter( ~/.*${tumorpattern}.*\.bam$/ ).collect() // collect these into a list format so you can use the .join operator to string together for -I input
    tumor_ch_bai = bai_dir.filter( ~/.*${tumorpattern}.*\.bai$/ )
    normal_ch = bam_dir.filter( ~/.*${normalpattern}.*\.bam$/ )
    normal_ch_bai = bai_dir.filter( ~/.*${normalpattern}.*\.bai$/ )

    // Normal sample id channel: take full filename and grab the sample id. Need this for the -normal arg in mutect2
    sample_id = normal_ch.map { filePath -> 
        def fileName = filePath.baseName // get file name without extensions, "TCGA-00-0000-00-00-000-0000-00_tumor_otherinfo_otherinfo.bam"
        def sampleName = fileName.split('_') // split the file name into an array by underscores [TCGA-00-0000-00-00-000-0000-00, tumor, otherinfo, otherinfo.bam]
        def listSample = sampleName as List // convert to a list to perform list operations
        def samplename = listSample[0]// grab the first element which is always the sample ID based on the RENAME process
        return samplename}
    sample_id_ch = sample_id.first() // convert to a value channel using .first()

    // Run mutect2
    MUTECT2(tumor_ch, tumor_ch_bai, normal_ch, normal_ch_bai, chrom_ch, sample_id_ch, params.mutect_idx, params.mutect_idx_fai, params.mutect_idx_dict)
    
    // Merge and prepare VCF
    BGZIP(MUTECT2.out.vcf) // concatenation requires bgzip'd files 
    vcfs_ch = BGZIP.out.vcf.collect() // collect all bgzip vcf outputs into a channel
    split_vcf_index = BGZIP.out.index.collect() // collect all bgzip index outputs into a channel
    // concatenate, normalize, and sort the VCF
    PREPAREVCF(vcfs_ch, split_vcf_index, sample_id_ch, params.mutect_idx, params.mutect_idx_fai, params.mutect_idx_dict)
    unfiltered_vcf = PREPAREVCF.out.normalized
    unfiltered_vcf_index = PREPAREVCF.out.index
    
    // Merge stats file
    stats = MUTECT2.out.stats
    stats_ch = stats.collect()
    MERGESTATS(stats_ch, sample_id_ch)
    filter_stats = MERGESTATS.out

    // Merge f1r2 read orientation files 
    f1r2files = MUTECT2.out.f1r2
    f1r2_ch = f1r2files.collect()
    LEARNORIENTATION(f1r2_ch, sample_id_ch)
    orientationmodel = LEARNORIENTATION.out

    // Filter mutect2 calls
    FILTERMUTECT(unfiltered_vcf, unfiltered_vcf_index, params.mutect_idx, params.mutect_idx_fai, params.mutect_idx_dict, filter_stats, orientationmodel, segment_table, contam_table, sample_id_ch)
    filter_vcf = FILTERMUTECT.out
    
    // filter for passing variants
    PASS(filter_vcf, sample_id_ch)

    // filter for variants above specific thresholds 

    // Annotate with snpEff
    ANNOTATE(PASS.out, sample_id_ch)
    
}
