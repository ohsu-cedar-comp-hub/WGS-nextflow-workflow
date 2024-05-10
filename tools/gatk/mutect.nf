#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Define the process for running MuTect2
process mutect2 {
    // Set maximum memory
    memory '40 GB'

    // Set output directory for MuTect2 results
    publishDir "${params.outdir}/svc", mode: 'copy'

    // Define input and output
    input:
    path tumor_bam_sorted
    path normal_bam_sorted
    path mutect_idx
    path pon
    path gnomad

    output:
    path "${tumor_bam.baseName}_unfiltered.vcf"
    path "${tumor_bam.baseName}_unfiltered.vcf.stats"
    path "${tumor_bam.baseName}_f1r2.tar.gz"

    // MuTect2 command
    script:
    """
    gatk Mutect2 \
        -R ${mutect_idx} \
        -I ${tumor_bam_sorted} \
        -I ${normal_bam_sorted} \
        --panel-of-normals ${pon} \
        --germline-resource ${gnomad} \
        -O ${tumor_bam_sorted.baseName}_unfiltered.vcf \
        --f1r2-tar-gz ${tumor_bam_sorted.baseName}_f1r2.tar.gz \
        -stats ${tumor_bam.baseName}_unfiltered.vcf.stats
    """
}