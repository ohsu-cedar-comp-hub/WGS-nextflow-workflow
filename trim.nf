#!/usr/bin/env nextflow

process trimmomatic {
    // Specify the directory where output files will be published
    publishDir "${params.outdir}", mode: 'copy'

    // Input: Raw FASTQ.gz file
    input:
    file(infile)

    // Output: Trimmed FASTQ.gz file
    output:
    file("${params.outdir}/${infile.baseName}_trimmed.fastq.gz")

    // Script to execute Trimmomatic
    script:
    """
    trimmomatic SE -phred33 ${infile} \
    ${params.outdir}/${infile.baseName}_trimmed.fastq.gz \
    ILLUMINACLIP:/home/exacloud/gscratch/CEDAR/grieco/nextflow/TruSeq3-SE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

workflow {
    // Run Trimmomatic for the input file
    trimmomatic(file(params.infile))
}