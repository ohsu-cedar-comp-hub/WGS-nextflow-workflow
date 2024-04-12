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
        -I ${params.tumor_bam} \
        -I ${params.normal_bam} \
        --panel-of-normals ${params.pon} \
        -O ${tumor_bam.baseName}_unfiltered.vcf \
        --f1r2-tar-gz ${tumor_bam.baseName}_f1r2.tar.gz
    """
}
