#!/usr/bin/env nextflow

// Process for getting pileup summaries
process GETPILEUPSUMMARIES {

    container "${params.container_gatk}"

    publishDir "${params.outdir}/summaries", mode: 'copy'

    input:
    path tumor_bam_sorted
    path normal_bam_sorted
    path exac

    output:
    path ("${tumor_bam_sorted.simpleName}.getpileupsummaries.table"), emit: tumor
    path("${normal_bam_sorted.simpleName}.getpileupsummaries.table"), emit: normal
    
    script:
    """
    gatk GetPileupSummaries \\
    -I ${tumor_bam_sorted} \\
    -V ${params.exac} \\
    -L ${params.exac} \\
    -O ${tumor_bam_sorted.simpleName}.getpileupsummaries.table

    gatk GetPileupSummaries \\
    -I ${normal_bam_sorted} \\
    -V ${params.exac} \\
    -L ${params.exac} \\
    -O ${normal_bam_sorted.simpleName}.getpileupsummaries.table
    """
}
