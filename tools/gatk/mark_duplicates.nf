#!/usr/bin/env nextflow

process MARKDUPLICATES {
    
    container "${params.container_gatk}"

    publishDir "${params.outdir}/metrics", mode: 'copy', pattern: "*_metrics.txt"
    
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

