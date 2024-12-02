#!/usr/bin/env nextflow

// Nextflow Pipeline Version
params.release = ""
params.releasedate = ""
params.githublink = "https://github.com/ohsu-cedar-comp-hub/WGS-nextflow-workflow/releases/tag/"

// Import tools
include { TUMORONLYGETPILEUPSUMMARIES } from '../../tools/gatk/get_pileup_summaries.nf'
include { TUMORONLYCALCULATECONTAMINATION } from '../../tools/gatk/calculate_contamination.nf'
// include { GETINTERVALS } from '../../tools/samtools/get_intervals.nf'
include { TUMORONLYMUTECT2 } from '../../tools/gatk/mutect2.nf'
include { BGZIP; PREPAREVCF } from '../../tools/bcftools/prepareVCFs.nf'
include { MERGESTATS } from '../../tools/bcftools/combineMutectStats.nf'
include { LEARNORIENTATION } from '../../tools/bcftools/combineF1R2files.nf'
include { FILTERMUTECT } from '../../tools/gatk/filter_mutect.nf'
include { REHEADER } from '../../tools/bcftools/reheader.nf'
include { FUNCOTATOR } from '../../tools/gatk/funcotator.nf'
include { SNPEFF } from '../../tools/snpeff/annotate_variants.nf'
include { PASS } from '../../tools/snpeff/sift_variants.nf'
include { ADDFILTER } from '../../tools/bcftools/filterVCF.nf'

tumor_ch = Channel.fromPath("${params.bam_files}/*.bam")
tumor_ch_bai = Channel.fromPath("${params.bam_files}/*.bai")

chromosomes = (1..22).collect { it.toString() } + ['X']
chrom_strings = Channel.from(chromosomes)
chrom_ch = chrom_strings.map { it -> "chr" + it }

// Begin main workflow
workflow {
    TUMORONLYMUTECT2(tumor_ch, tumor_ch_bai, chrom_ch, params.sample_id, params.mutect_idx, params.mutect_idx_fai, params.mutect_idx_dict, params.pon_vcf, params.pon_tbi, params.pon_idx, params.pon_tar)
    
    // gatk getpileupsummaries
    TUMORONLYGETPILEUPSUMMARIES(tumor_ch, tumor_ch_bai, params.exac)
    tumor_table = TUMORONLYGETPILEUPSUMMARIES.out

    // gatk calculate contamination from pileup summaries
    TUMORONLYCALCULATECONTAMINATION(tumor_table)
    contam_table = TUMORONLYCALCULATECONTAMINATION.out.contamination.collect()
    segment_table = TUMORONLYCALCULATECONTAMINATION.out.segment.collect()

    // Merge and prepare VCF
    BGZIP(TUMORONLYMUTECT2.out.vcf) // concatenation requires bgzip'd files 
    vcfs_ch = BGZIP.out.vcf.collect() // collect all bgzip vcf outputs into a channel
    split_vcf_index = BGZIP.out.index.collect() // collect all bgzip index outputs into a channel
    // concatenate, normalize, and sort the VCF
    PREPAREVCF(vcfs_ch, split_vcf_index, params.sample_id, params.mutect_idx, params.mutect_idx_fai, params.mutect_idx_dict)
    unfiltered_vcf = PREPAREVCF.out.vcf
    unfiltered_vcf_index = PREPAREVCF.out.index
    
    // Merge stats file
    stats = TUMORONLYMUTECT2.out.stats
    stats_ch = stats.collect()
    MERGESTATS(stats_ch, params.sample_id)
    filter_stats = MERGESTATS.out

    // Merge f1r2 read orientation files 
    f1r2files = TUMORONLYMUTECT2.out.f1r2
    f1r2_ch = f1r2files.collect()
    LEARNORIENTATION(f1r2_ch, params.sample_id)
    orientationmodel = LEARNORIENTATION.out

    // Filter mutect2 calls
    FILTERMUTECT(unfiltered_vcf, unfiltered_vcf_index, params.mutect_idx, params.mutect_idx_fai, params.mutect_idx_dict, filter_stats, orientationmodel, segment_table, contam_table, params.sample_id)
    
    // Add nextflow workflow versioning to VCF header
    REHEADER(FILTERMUTECT.out)

    // filter for passing variants
    PASS(REHEADER.out, params.sample_id)

    // filter for variants above certain allelic depth, VAF, etc using bcftools
    ADDFILTER(PASS.out, params.sample_id)
    vcf = ADDFILTER.out
    
    // Annotate with funcotator
    FUNCOTATOR(vcf, 
        params.mutect_idx, params.mutect_idx_fai, params.mutect_idx_dict,
        params.funcotator_data,
        params.sample_id)

    // Annotate with snpEff
    SNPEFF(vcf, params.sample_id)
}
