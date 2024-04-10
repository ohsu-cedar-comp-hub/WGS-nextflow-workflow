#!/usr/bin/env nextflow

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
