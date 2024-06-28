#!/usr/bin/env nextflow

process ASCAT {
    container "${params.container_cgpwgs}"
    publishDir "${params.outdir}/cgp-wgs/ascat", mode: 'copy'

    input:
    tumor bam
    normal bam

    output:
    tumor copy number bed
    normal copy number bed

    script:
    """
    """
}