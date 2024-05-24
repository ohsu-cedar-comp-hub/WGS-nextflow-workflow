#!/usr/bin/env nextflow

process MULTIQC {

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    // require all files from fastqc before running multiqc 
    file("*")    

    output:
    file("multiqc_report.html")

    script:
    """
    multiqc .
    """
}
