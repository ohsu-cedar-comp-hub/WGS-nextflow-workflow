#!/usr/bin/env nextflow

process SortMarkedDuplicates {
    // Set maximum memory
    memory '40 GB'
    
    publishDir "${params.outdir}/aligned/markduplicates/sorted", mode: 'copy'

    input:
    path bam_duplicates_unsorted
    val ID

    output: 
    file("${bam_duplicates_unsorted.baseName}_sorted.bam")
    file("${bam_duplicates_unsorted.baseName}_sorted.bam.bai")


    script:
    """
    samtools sort ${bam_duplicates_unsorted} > ${bam_duplicates_unsorted.baseName}_sorted.bam
    samtools index ${bam_duplicates_unsorted.baseName}_sorted.bam > ${bam_duplicates_unsorted.baseName}_sorted.bam.bai
    """
}