#!/usr/bin/env nextflow

process TOFASTQSORT {
    
    publishDir "${params.outdir}/tofastq"

    container "${params.container_samtools}"

    input:
    path bam_cram

    output:
    path("*.fastq.gz")

    script:
    """
    touch samtools_test.bam
    """
}

process TOFASTQ {
        
    publishDir "${params.outdir}/tofastq"

    container "${params.container_samtools}"

    input:
    path bam_cram

    output:
    path("*.fastq.gz")

    script:

    """
    samtools collate -u -O ${bam_cram} | samtools fastq -1 ${bam_cram.baseName}_R1.fastq.gz -2 ${bam_cram.baseName}_R2.fastq.gz -0 /dev/null -s /dev/null -n
    """
}