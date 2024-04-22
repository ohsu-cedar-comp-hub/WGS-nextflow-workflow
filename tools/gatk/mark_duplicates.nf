#!/usr/bin/env nextflow

process MarkDuplicates {
    publishDir "${params.outdir}/aligned/markduplicates", mode: 'copy'

    input:
    path tumor_bam_sorted
    path normal_bam_sorted
    val ID

    output: 
    file "${tumor_bam_sorted.baseName}_marked_duplicates.bam"
    file "${tumor_bam_sorted.baseName}_marked_duplicates_metrics.txt"
    file "${normal_bam_sorted.baseName}_marked_duplicates.bam"
    file "${normal_bam_sorted.baseName}_marked_duplicates_metrics.txt"

    script:
    """
    gatk MarkDuplicates -I ${tumor_bam_sorted} \
        -O ${tumor_bam_sorted.baseName}_marked_duplicates.bam \
        -M ${tumor_bam_sorted.baseName}_marked_duplicates_metrics.txt \
        --VALIDATION_STRINGENCY LENIENT

   gatk MarkDuplicates -I ${normal_bam_sorted} \
        -O ${normal_bam_sorted.baseName}_marked_duplicates.bam \
        -M ${normal_bam_sorted.baseName}_marked_duplicates_metrics.txt \
        --VALIDATION_STRINGENCY LENIENT
    """

}

