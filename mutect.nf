#!/usr/bin/env nextflow


// Define the process for running MuTect2
process mutect2 {
    // Set output directory for MuTect2 results
    publishDir "${params.outdir}/svc", mode: 'copy'

    // Define input and output
    input:
    path bam_sorted
    path mutect_idx
    path pon

    output:
    file("${bam_sorted.baseName}_svc.vcf")

    // MuTect2 command
    script:
    """
    gatk Mutect2 \
        -R ${params.mutect_idx} \
        -I ${params.bam_sorted} \
        --panel-of-normals ${params.pon} \
        -O ${bam_sorted.baseName}_svc.vcf
    """
}

// Define the workflow
workflow {
    // define input paramaters for Mutect2
    input_file=file(params.bam_sorted)
    idx=file(params.mutect_idx)
    pon=file(params.pon)
    // Run Mutect2
    mutect2(input_file, idx, pon)
}