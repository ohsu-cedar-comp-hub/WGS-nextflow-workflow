#!/usr/bin/env nextflow

process TOFASTQ {
    
    publishDir "${params.outdir}/preptofastq"

    container "${params.container_samtools}"

    input:
    path bam_cram

    output:
    path("*.fastq.gz")

/*
    script:
    """
    samtools sort ${bam_cram} > ${bam_cram.baseName}_sorted.bam
    samtools collate -u -O ${bam_cram.baseName}_sorted.bam | samtools fastq -1 ${bam_cram.baseName}_R1.fastq.gz -2 ${bam_cram.baseName}_R2.fastq.gz
    """
*/
    script:
    """
    touch ${bam_cram.baseName}_sorted.bam
    touch ${bam_cram.baseName}_R1.fastq.gz
    touch ${bam_cram.baseName}_R2.fastq.gz
    """
}

