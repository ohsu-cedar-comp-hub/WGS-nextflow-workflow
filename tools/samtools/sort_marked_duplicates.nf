#!/usr/bin/env nextflow

process SortMarkedDuplicates {
    publishDir "${params.outdir}/aligned/markduplicates/sorted", mode: 'copy'

    input:
    path bam_duplicates_unsorted

    output: 
    file("${bam_duplicates_unsorted.baseName}_sorted.bam")


    script:
    """
    samtools sort ${bam_duplicates_unsorted} > ${bam_duplicates_sorted.baseName}_sorted.bam
    """
}