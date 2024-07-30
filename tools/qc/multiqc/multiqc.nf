#!/usr/bin/env nextflow

process MULTIQC {

    container "${params.container_multiqc}"

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    // require all files from fastqc before running multiqc 
    file("*")
    tuple val(sample_id), path(reads)

    samplename = sample_id.first()

    output:
    file("${samplename}_multiqc_report.html")

    script:
    """
    multiqc --filename ${samplename}_multiqc_report.html .
    """
}
