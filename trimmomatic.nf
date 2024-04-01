#!/usr/bin/env nextflow

process trimmomaticPE {
    publishDir "${params.outdir}/trim_reads", mode: 'copy'

    input:
    path read1
    path read2
    path truseq3pefile

    output:
    file("${read1.baseName}_1P.fastq.gz")
    file("${read1.baseName}_1U.fastq.gz")
    file("${read2.baseName}_2P.fastq.gz")
    file("${read2.baseName}_2U.fastq.gz")

    script:
    """
     trimmomatic \
     PE -phred33 \
     ${read1} \
     ${read2} \
     ${read1.baseName}_1P.fastq.gz \
     ${read1.baseName}_1U.fastq.gz \
     ${read2.baseName}_2P.fastq.gz \
     ${read2.baseName}_2U.fastq.gz \
     ILLUMINACLIP:"${truseq3pefile}":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
    """
}

workflow {
    read1=file(params.read1)
    read2=file(params.read2)
    truseq3pefile=file(params.truseq3pefile)

    trimmomaticPE(read1, read2, truseq3pefile)
}
