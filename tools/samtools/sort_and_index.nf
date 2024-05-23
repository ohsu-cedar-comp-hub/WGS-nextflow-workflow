#!/usr/bin/env nextflow

process SORT {
    // Set maximum memory
    memory '40 GB'
    
    publishDir "${params.outdir}/aligned/sort_index", mode: 'copy'
    
    input:
    path bam_unsorted

    output: 
    file("${bam_unsorted.baseName}_sorted.bam")
    
    script:
    """
    samtools sort ${bam_unsorted} > ${bam_unsorted.baseName}_sorted.bam
    """
}

process SORTANDINDEX {
    // Set maximum memory
    memory '40 GB'
    
    publishDir "${params.outdir}/aligned/markduplicates/sorted", mode: 'copy'
    
    input:
    path bam_unsorted

    output: 
    path("${bam_unsorted.baseName}_sorted.bam"), emit: bam
    path("${bam_unsorted.baseName}_sorted.bam.bai"), emit: bai
    
    script:
    """
    samtools sort ${bam_unsorted} > ${bam_unsorted.baseName}_sorted.bam
    samtools index ${bam_unsorted.baseName}_sorted.bam > ${bam_unsorted.baseName}_sorted.bam.bai
    """
}