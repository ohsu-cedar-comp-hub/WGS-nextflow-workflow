#!/usr/bin/env nextflow

// Define the process for running MuTect2
process MUTECT2 {
    // Set maximum memory
    //memory '40 GB'

    cpus 1 // set cpu to 1: gatk discourages multithreading

    container "${params.container_gatk}"

    input:
    path tumor_bam_sorted 
    path normal_bam_sorted
    val chrom 
    val sample_id

    output:
    path "${sample_id}_${chrom}_unfiltered.vcf", emit: vcf
    path "${sample_id}_${chrom}_f1r2.tar.gz", emit: f1r2
    path "${sample_id}_${chrom}_unfiltered.vcf.stats", emit: stats
    path "${sample_id}_${chrom}_unfiltered.vcf.idx", emit: index

    script:
    """
    gatk Mutect2 \
    -R ${params.mutect_idx} \
        -I ${tumor_bam_sorted} \
        -I ${normal_bam_sorted} \
        --panel-of-normals ${params.pon} \
        -L ${chrom} \
        --germline-resource ${params.gnomad} \
        -O ${sample_id}_${chrom}_unfiltered.vcf \
        --f1r2-tar-gz ${sample_id}_${chrom}_f1r2.tar.gz
    """
}
