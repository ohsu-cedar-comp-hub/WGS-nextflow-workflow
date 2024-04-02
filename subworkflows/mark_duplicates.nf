#!/usr/bin/env nextflow

process MarkDuplicates {
    publishDir "${params.outdir}/aligned/markduplicates", mode: 'copy'

    input:
    path bam_sorted

    output: 
    file "${bam_sorted.baseName}_marked_duplicates.bam"
    file "${bam_sorted.baseName}_marked_duplicates_metrics.txt"

    script:
    """
    gatk MarkDuplicates \
        -I ${bam_sorted} \
        -O ${bam_sorted.baseName}_marked_duplicates.bam \
        -M ${bam_sorted.baseName}_marked_duplicates_metrics.txt \
        --VALIDATION_STRINGENCY LENIENT
    """

}
