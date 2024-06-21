#!/usr/bin/env nextflow

process MULTIQC {

    container "${params.container_multiqc}"

    publishDir "${params.outdir}/multiqc"

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
