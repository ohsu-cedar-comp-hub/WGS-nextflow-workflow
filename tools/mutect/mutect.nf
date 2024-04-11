#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Define the process for running MuTect2
process mutect2 {
    // Set output directory for MuTect2 results
    publishDir "${params.outdir}/svc", mode: 'copy'

    // Define input and output
    input:
    path tumor_bam
    path normal_bam
    path mutect_idx
    path pon

    output:
    path "${tumor_bam.baseName}_unfiltered.vcf"
    path "${tumor_bam.baseName}_unfiltered.vcf.stats"
    path "${tumor_bam.baseName}_f1r2.tar.gz"

    // MuTect2 command
    script:
    """
    gatk Mutect2 \
        -R ${params.mutect_idx} \
        -I ${params.bam_sorted} \
        --panel-of-normals ${params.pon} \
        -O ${tumor_bam.baseName}_unfiltered.vcf \
        --f1r2-tar-gz ${tumor_bam.baseName}_f1r2.tar.gz
    """
}

// Define the workflow
workflow {
    // define input paramaters for Mutect2
    tumor_bam=file(params.tumor_bam)
    normal_bam=file(params.normal_bam)
    idx=file(params.mutect_idx)
    pon=file(params.pon)
    // Run Mutect2
    mutect2(input_file, idx, pon)
}