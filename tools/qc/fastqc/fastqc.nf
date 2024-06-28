#!/usr/bin/env nextflow

// Define the process for running FastQC

process FASTQC {

    container "${params.container_fastqc}"

    input:
    tuple val(sample_id), path(reads)
    path outdir

    output:
    path("${sample_id}*_fastqc.zip"), emit: zip
    path("${sample_id}*_fastqc.html"), emit: html
    

    script:
    """
    # run fastqc on pair of files
    fastqc ${reads[0]} ${reads[1]}
    """
}
