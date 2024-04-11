#!/usr/bin/env nextflow

process sortAndIndex {
    publishDir "${params.outdir}/aligned/sort_index", mode: 'copy'

    input:
    path bam_unsorted

    output: 
    file("${bam_unsorted.baseName}_sorted_indexed.bam")
    file("${bam_unsorted.baseName}_sorted_indexed.bam.bai")

    script:
    """
    samtools sort ${bam_unsorted} > ${bam_unsorted.baseName}_sorted_indexed.bam
    samtools index ${bam_unsorted.baseName}_sorted_indexed.bam > ${bam_unsorted.baseName}_sorted_indexed.bam.bai
    """
}
