#!/usr/bin/env nextflow

process TRIMMOMATICPE {
    // Set maximum memory
    // memory '40 GB'

    publishDir "${params.outdir}/trim_reads", mode: 'copy'

    input:
    tuple val(sample_id), val(reads)
    path truseq3pefile
    path outdir

    output:
    tuple val(sample_id), path("*P.fastq.gz"), emit: trim_reads
    path("*U.fastq.gz"), optional:true, emit: unpaired_reads

    script:
    """
    java -jar /bin/trimmomatic.jar \
     PE -phred33 \
     ${reads[0]} \
     ${reads[1]} \
     ${reads[0].simpleName}_1P.fastq.gz \
     ${reads[0].simpleName}_1U.fastq.gz \
     ${reads[1].simpleName}_2P.fastq.gz \
     ${reads[1].simpleName}_2U.fastq.gz \
     ILLUMINACLIP:"${truseq3pefile}":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
    """
}