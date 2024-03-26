#!/usr/bin/env nextflow

// Define the process for trimmomatic
process trimmomaticPE {
    publishDir "${params.outdir}", mode: 'copy'

     // Define input and output
    input:
    path Read1
    path Read2
    path truSeq3PeFile

    output:
    file("${Read1.baseName}_1P.fastq.gz")
    file("${Read1.baseName}_1U.fastq.gz")
    file("${Read2.baseName}_2P.fastq.gz")
    file("${Read2.baseName}_2U.fastq.gz")

    // trimmomatic command
    script:
    """
     trimmomatic \
     PE -phred33 \
     ${Read1} \
     ${Read2} \
     ${Read1.baseName}_1P.fastq.gz \
     ${Read1.baseName}_1U.fastq.gz \
     ${Read2.baseName}_2P.fastq.gz \
     ${Read2.baseName}_2U.fastq.gz \
     ILLUMINACLIP:"${truSeq3PeFile}":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
    """
}
