#!/usr/bin/env nextflow

process sortAndIndex {
    // Set maximum memory
    memory '40 GB'

    publishDir "${params.outdir}/aligned/sort_index", mode: 'copy'

    input:
    path bam_unsorted
    val ID

    output: 
    file("${bam_unsorted.baseName}_sorted_indexed.bam")
    file("${bam_unsorted.baseName}_sorted_indexed.bam.bai")

    script:
    """
    samtools sort ${bam_unsorted} > ${bam_unsorted.baseName}_sorted_indexed.bam
    samtools index ${bam_unsorted.baseName}_sorted_indexed.bam > ${bam_unsorted.baseName}_sorted_indexed.bam.bai
    """
}
