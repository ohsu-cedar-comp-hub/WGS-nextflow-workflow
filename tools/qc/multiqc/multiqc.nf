#!/usr/bin/env nextflow

process MULTIQC {

    container "${params.container_multiqc}"

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    // require all files from fastqc before running multiqc 
    file("*")
    tuple val(sample_id), val(type)

    output:
    file("${sample_id}_multiqc_report.html")

    script:
    """
    multiqc --filename ${sample_id}_multiqc_report.html .
    """
}
