#!/usr/bin/env nextflow

// Define the process for running FastQC

process FASTQC {
    // debug true

    conda "${params.conda_env}"
    // container "${params.container_fastqc}"

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
    fastqc ${reads[0]} ${reads[1]}
    """
}

all_pairs_ch = Channel.fromFilePairs(params.all_read_pairs)

workflow {
    FASTQC(all_pairs_ch, params.outdir)
}