#!/usr/bin/env nextflow

// Define the process for running MuTect2
process MUTECT2 {
    // Set maximum memory
    // memory '40 GB'

    cpus 1 // set cpu to 1: gatk discourages multithreading

    container "${params.container_gatk}"

    input:
    path tumor_bam_sorted 
    path tumor_bam_sorted_bai
    path normal_bam_sorted
    path normal_bam_sorted_bai
    val chrom 
    val sample_id
    path mutect_idx
    path mutect_idx_fai
    path mutect_idx_dict

    output:
    path "${sample_id}_${chrom}_unfiltered.vcf", emit: vcf
    path "${sample_id}_${chrom}_f1r2.tar.gz", emit: f1r2
    path "${sample_id}_${chrom}_unfiltered.vcf.stats", emit: stats
    path "${sample_id}_${chrom}_unfiltered.vcf.idx", emit: index

    script:

    // normal name is whatever the SM: label in your @RG is

    """
    gatk Mutect2 \
    -R ${mutect_idx} \
        -I ${tumor_bam_sorted.join(' -I ')} \
        -I ${normal_bam_sorted} \
        -normal ${sample_id} \
        --panel-of-normals ${params.pon} \
        -L ${chrom} \
        --germline-resource ${params.gnomad} \
        -O ${sample_id}_${chrom}_unfiltered.vcf \
        --f1r2-tar-gz ${sample_id}_${chrom}_f1r2.tar.gz
    """
}
