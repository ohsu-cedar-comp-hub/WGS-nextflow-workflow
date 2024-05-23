#!/usr/bin/env nextflow

// Define the process for running FastQC

process FASTQC {
    // debug true
    
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    path reads
    path outdir

    output:
    path("*.zip"), emit: zip
    path("*.html"), emit: html
    
    script:
    """
    # fastqc requires a pre-made directory. check it exists. if not, make it.
    if [ -d $outdir/fastqc ]; then
        :
    else
        mkdir $outdir/fastqc
    fi 
    
    # run fastqc
    /usr/local/FastQC/fastqc $reads
    """
}

