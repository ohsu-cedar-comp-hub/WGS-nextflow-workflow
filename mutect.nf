#!/usr/bin/env nextflow


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
    file("${tumor_bam.baseName}.vcf")

    // MuTect2 command
    script:
    """
    gatk Mutect2 \
        -R ${params.mutect_idx} \
        -I ${params.tumor_bam} \
        -I ${params.normal_bam} \
        --f1r2-tar-gz f1r2.tar.gz \
        --panel-of-normals ${params.pon} \
        -O ${bam_sorted.baseName}_svc.vcf
    """
}

// Define the workflow
workflow {
    // define input paramaters for Mutect2
    input_normal_file=file(params.normal_bam)
    input_tumor_file=file(params.tumor_bam)
    idx=file(params.mutect_idx)
    pon=file(params.pon)
    // Run Mutect2
    mutect2(input_normal_file, input_tumor_file, idx, pon)
}
