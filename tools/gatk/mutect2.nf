#!/usr/bin/env nextflow

process MUTECT2 {

    maxForks 12 // set this when running on local scratch to parallelize; up to the max number of cpus available minus 1
    cpus 1 // set cpu to 1: gatk discourages multithreading
    container "${params.container_gatk}"

    input:
    path tumor_bam
    path tumor_bam_sorted_bai // only necessary when nextflow can't resolve path from symlink
    path normal_bam
    path normal_bam_sorted_bai // only necessary when nextflow can't resolve path from symlink
    each chrom // repeat this process for each item in the chrom channel
    val sample_id 
    val normal_command
    path mutect_idx // nextflow can't resolve the rest of these files from symlink to mutect_idx, input paths to each
    path mutect_idx_fai
    path mutect_idx_dict
    path pon_vcf // nextflow can't resolve the rest of these files from symlink to pon_vcf, input paths to each
    path pon_tbi
    path pon_idx
    path pon_tar

    output:
    path "${sample_id}_${chrom}_unfiltered.vcf", emit: vcf
    path "${sample_id}_${chrom}_f1r2.tar.gz", emit: f1r2
    path "${sample_id}_${chrom}_unfiltered.vcf.stats", emit: stats
    path "${sample_id}_${chrom}_unfiltered.vcf.idx", emit: index

    script:

    // normal name (sample_id) is whatever the SM: label in your @RG is in the bam file

    """
    gatk Mutect2 \
    -R ${mutect_idx} \
        -I ${tumor_bam.join(' -I ')}  \
        -I ${normal_bam.join(' -I ')} \
        ${normal_command} \
        --panel-of-normals ${params.pon_vcf} \
        -L ${chrom} \
        --germline-resource ${params.gnomad} \
        -O ${sample_id}_${chrom}_unfiltered.vcf \
        --f1r2-tar-gz ${sample_id}_${chrom}_f1r2.tar.gz
    """
}

// Define the process for running MuTect2
process TUMORONLYMUTECT2 {
    
    maxForks 12 // set this when running on local scratch to parallelize; up to the max number of cpus available minus 1
    cpus 1 // set cpu to 1: gatk discourages multithreading
    container "${params.container_gatk}"

    input:
    path tumor_bam
    path tumor_bam_sorted_bai // only necessary when nextflow can't resolve path from symlink
    each chrom // repeat this process for each item in the chrom channel
    val sample_id 
    path mutect_idx // nextflow can't resolve the rest of these files from symlink to mutect_idx, input paths to each
    path mutect_idx_fai
    path mutect_idx_dict
    path pon_vcf // nextflow can't resolve the rest of these files from symlink to pon_vcf, input paths to each
    path pon_tbi
    path pon_idx
    path pon_tar

    output:
    path "${sample_id}_${chrom}_unfiltered.vcf", emit: vcf
    path "${sample_id}_${chrom}_f1r2.tar.gz", emit: f1r2
    path "${sample_id}_${chrom}_unfiltered.vcf.stats", emit: stats
    path "${sample_id}_${chrom}_unfiltered.vcf.idx", emit: index

    script:

    """
    gatk Mutect2 -R ${mutect_idx} -I ${tumor_bam} --panel-of-normals ${params.pon_vcf} -L ${chrom} --germline-resource ${params.gnomad} -O ${sample_id}_${chrom}_unfiltered.vcf --f1r2-tar-gz ${sample_id}_${chrom}_f1r2.tar.gz
    """
}
