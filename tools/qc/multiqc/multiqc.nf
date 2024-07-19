#!/usr/bin/env nextflow

process MULTIQC {

    container "${params.container_multiqc}"

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    // require all files from fastqc before running multiqc 
    file("*")
    tuple val(sampleid), path(reads)    

    output:
    file("${reads[0].basename}_multiqc_report.html")

    script:
    """
    multiqc --filename ${reads[0].basename}_multiqc_report.html .
    """
}
