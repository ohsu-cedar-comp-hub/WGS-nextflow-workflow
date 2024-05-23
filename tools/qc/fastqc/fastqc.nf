#!/usr/bin/env nextflow

// Define the process for running FastQC

process FASTQC {
    // debug true

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path outdir

    output:
    path("${sample_id}*_fastqc.zip"), emit: zip
    path("${sample_id}*_fastqc.html"), emit: html
    

    script:
    """
    # fastqc requires a pre-made directory. check it exists. if not, make it.
    if [ -d $outdir/fastqc ]; then
        :
    else
        mkdir $outdir/fastqc
    fi 
    
    # run fastqc on pair of files
    /usr/local/FastQC/fastqc ${reads[0]} ${reads[1]}
    """
}

