#!/usr/bin/env nextflow

process SORT {
    // Set maximum memory
    // memory '40 GB'
    
    container "${params.container_samtools}"

    publishDir "${params.outdir}/aligned/sort", mode: 'copy'
    
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
<<<<<<< HEAD
    // memory "${params.memory}"
=======
    // memory '40 GB'
>>>>>>> remove memory requirements for testing
    
    container "${params.container_samtools}"

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