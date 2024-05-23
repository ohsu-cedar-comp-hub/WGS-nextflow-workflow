#!/usr/bin/env nextflow

//Define process for multiQC
process MULTIQC {
    debug true 

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path read

    output:
    path("${read.simpleName}_multiqc_report.html")

    script:
    """
    multiqc $read -o ${read.simpleName}_multiqc_report.html
    """
}
