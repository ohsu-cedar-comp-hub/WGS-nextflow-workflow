#!/usr/bin/env nextflow

process trimmomaticPE {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path Read1
    path Read2

    output:
    file "${params.outdir}/${Read1.baseName}_1P.fastq.gz"
    file "${params.outdir}/${Read1.baseName}_1U.fastq.gz"
    file "${params.outdir}/${Read2.baseName}_2P.fastq.gz"
    file "${params.outdir}/${Read2.baseName}_2U.fastq.gz"

    script:
    """
    trimmomatic PE -phred33 \
        ${Read1} ${Read2} \
        ${params.outdir}/${Read1.baseName}_1P.fastq.gz ${params.outdir}/${Read1.baseName}_1U.fastq.gz \
        ${params.outdir}/${Read2.baseName}_2P.fastq.gz ${params.outdir}/${Read2.baseName}_2U.fastq.gz \
        ILLUMINACLIP:"/home/exacloud/gscratch/CEDAR/grieco/nextflow/TruSeq3-PE.fa":2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

workflow {
    Read1=file(params.Read1)
    Read2=file(params.Read2)

    trimmomaticPE(Read1, Read2)
}