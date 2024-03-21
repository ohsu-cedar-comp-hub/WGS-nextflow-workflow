#!/usr/bin/env nextflow

process trimmomaticPE {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path Read1
    path Read2
    path truSeq3PeFile

    output:
    file("${Read1.baseName}_1P.fastq.gz")
    file("${Read1.baseName}_1U.fastq.gz")
    file("${Read2.baseName}_2P.fastq.gz")
    file("${Read2.baseName}_2U.fastq.gz")

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

workflow {
    Read1=file(params.Read1)
    Read2=file(params.Read2)
    truSeq3PeFile=file(params.truSeq3PeFile)

    trimmomaticPE(Read1, Read2, truSeq3PeFile)
}
