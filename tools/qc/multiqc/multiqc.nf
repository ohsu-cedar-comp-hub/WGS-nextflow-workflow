#!/usr/bin/env nextflow

process MULTIQC {

    container "${params.container_multiqc}"

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    // require all files from fastqc before running multiqc 
    file("*")
    val sampleid  

    output:
    file("${sampleid}_multiqc_report.html")

    script:
    """
    multiqc --filename ${sampleid}_multiqc_report.html .
    """
}
