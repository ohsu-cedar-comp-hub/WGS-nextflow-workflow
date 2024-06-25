#!/usr/bin/env nextflow

process TOFASTQSORT {
    
    publishDir "${params.outdir}/tofastq", mode: 'copy'

    container "${params.container_samtools}"

    input:
    path bam_cram

    output:
    path("*.${params.filesuffix}")

    script:
    """
    samtools sort ${bam_cram} > ${bam_cram.baseName}_sorted.${params.filesuffix}
    """
}

process TOFASTQ {
        
    publishDir "${params.outdir}/tofastq", mode: 'copy'

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