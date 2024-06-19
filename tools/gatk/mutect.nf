#!/usr/bin/env nextflow

// Define the process for running MuTect2
process MUTECT2 {
    // Set maximum memory
    //memory '40 GB'

    // Set output directory for MuTect2 results
    // publishDir "${params.outdir}/svc", mode: 'copy'  // dont include this only copy the final file to save space. the split vcfs are still available in the work directory.

    container "${params.container_gatk}"

    // Define input and output
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

    // MuTect2 command
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
