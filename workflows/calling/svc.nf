#!/usr/bin/env nextflow

// Nextflow Pipeline Version
params.release = "v0.2.4"
params.releasedate = ""
params.githublink = "https://github.com/ohsu-cedar-comp-hub/WGS-nextflow-workflow/releases/tag/v0.2.4"

// Create queue channels for BAM and BAI files (consumable) from the BAM file directory
bam_dir = Channel.fromPath("${params.bam_files}/*.bam")
bai_dir = Channel.fromPath("${params.bam_files}/*.bai")

// Define the list of chromosomes + create a channel emitting each chromosome. Pass this channel instead of an intervals file for splitting Mutect2 jobs

// Import tools
include { GETPILEUPSUMMARIES } from '../../tools/gatk/get_pileup_summaries.nf'
include { CALCULATECONTAMINATION } from '../../tools/gatk/calculate_contamination.nf'
include { GETINTERVALS } from '../../tools/samtools/get_intervals.nf'
include { MUTECT2 } from '../../tools/gatk/mutect2.nf'
include { BGZIP; PREPAREVCF } from '../../tools/bcftools/prepareVCFs.nf'
include { MERGESTATS } from '../../tools/bcftools/combineMutectStats.nf'
include { LEARNORIENTATION } from '../../tools/bcftools/combineF1R2files.nf'
include { FILTERMUTECT } from '../../tools/gatk/filter_mutect.nf'
include { REHEADER } from '../../tools/bcftools/reheader.nf'
include { FUNCOTATOR } from '../../tools/gatk/funcotator.nf'
include { SNPEFF } from '../../tools/snpeff/annotate_variants.nf'
include { PASS } from '../../tools/snpeff/sift_variants.nf'
include { ADDFILTER } from '../../tools/bcftools/filterVCF.nf'

// Begin main workflow
workflow {
    // separate out tumor and normal samples into two different channels and collect into a list, so multiple files can be passed to Mutect2 in the case of several tumors/normals
    def tumorpattern = params.tumor_namepattern
    def normalpattern = params.normal_namepattern
    tumor_ch = bam_dir.filter( ~/.*${tumorpattern}.*\.bam$/ ).collect()
    tumor_ch_bai = bai_dir.filter( ~/.*${tumorpattern}.*\.bai$/ ).collect()
    normal_ch = bam_dir.filter( ~/.*${normalpattern}.*\.bam$/ ).collect()
    normal_ch_bai = bai_dir.filter( ~/.*${normalpattern}.*\.bai$/ ).collect()
    first_tumor = tumor_ch.flatten().first()
    first_tumor_bai = tumor_ch_bai.flatten().first()

    GETINTERVALS(first_tumor, first_tumor_bai)
    chrom_ch = GETINTERVALS.out.map { it.toString().split('\n') }.flatten().collect()
    chrom_ch.view()
    // Run Mutect2 as matched tumor(s)-normal(s)
    MUTECT2(tumor_ch, tumor_ch_bai,
        normal_ch, normal_ch_bai,
        chrom_ch, 
        params.normalsample_id, 
        params.normal_command,
        params.mutect_idx, params.mutect_idx_fai, params.mutect_idx_dict,
        params.pon_vcf, params.pon_tbi, params.pon_idx, params.pon_tar)
    
    // gatk getpileupsummaries
    GETPILEUPSUMMARIES(bam_dir, bai_dir, params.exac)
    tumor_table = GETPILEUPSUMMARIES.out.tumor
    normal_table = GETPILEUPSUMMARIES.out.normal.first() // assuming only one normal is used

    // gatk calculate contamination from pileup summaries
    CALCULATECONTAMINATION(tumor_table, normal_table)
    contam_table = CALCULATECONTAMINATION.out.contamination.collect()
    segment_table = CALCULATECONTAMINATION.out.segment.collect()

    // Merge and prepare VCF
    BGZIP(MUTECT2.out.vcf) // concatenation requires bgzip'd files 
    vcfs_ch = BGZIP.out.vcf.collect() // collect all bgzip vcf outputs into a channel
    split_vcf_index = BGZIP.out.index.collect() // collect all bgzip index outputs into a channel
    // concatenate, normalize, and sort the VCF
    PREPAREVCF(vcfs_ch, split_vcf_index, params.normalsample_id, params.mutect_idx, params.mutect_idx_fai, params.mutect_idx_dict)
    unfiltered_vcf = PREPAREVCF.out.vcf
    unfiltered_vcf_index = PREPAREVCF.out.index
    
    // Merge stats file
    stats = MUTECT2.out.stats
    stats_ch = stats.collect()
    MERGESTATS(stats_ch, params.normalsample_id)
    filter_stats = MERGESTATS.out

    // Merge f1r2 read orientation files 
    f1r2files = MUTECT2.out.f1r2
    f1r2_ch = f1r2files.collect()
    LEARNORIENTATION(f1r2_ch, params.normalsample_id)
    orientationmodel = LEARNORIENTATION.out

    // Filter mutect2 calls
    FILTERMUTECT(unfiltered_vcf, unfiltered_vcf_index, params.mutect_idx, params.mutect_idx_fai, params.mutect_idx_dict, filter_stats, orientationmodel, segment_table, contam_table, params.normalsample_id)
    
    // Add nextflow workflow versioning to VCF header
    REHEADER(FILTERMUTECT.out)

    // filter for passing variants
    PASS(REHEADER.out, params.normalsample_id)

    // filter for variants above certain allelic depth, VAF, etc using bcftools
    ADDFILTER(PASS.out, params.normalsample_id)
    vcf = ADDFILTER.out
    
    // Annotate with funcotator
    FUNCOTATOR(vcf, 
        params.mutect_idx, params.mutect_idx_fai, params.mutect_idx_dict,
        params.funcotator_data,
        params.normalsample_id)

    // Annotate with snpEff
    SNPEFF(vcf, params.normalsample_id)
}
