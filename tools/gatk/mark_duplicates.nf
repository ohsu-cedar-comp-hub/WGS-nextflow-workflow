#!/usr/bin/env nextflow

process MARKDUPLICATES {
    publishDir "${params.outdir}/aligned/markduplicates", mode: 'copy'

    input:
    path bam_sorted

    output: 
    path("${bam_sorted.baseName}_marked_duplicates.bam"), emit: bam
    path("${bam_sorted.baseName}_marked_duplicates_metrics.txt"), emit: metrics

    script:
    """
    gatk MarkDuplicates -I ${bam_sorted} \
        -O ${bam_sorted.baseName}_marked_duplicates.bam \
        -M ${bam_sorted.baseName}_marked_duplicates_metrics.txt \
        --VALIDATION_STRINGENCY LENIENT
    """

}

