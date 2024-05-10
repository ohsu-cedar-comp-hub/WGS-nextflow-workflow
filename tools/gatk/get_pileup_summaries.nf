#!/usr/bin/env nextflow

// Process for getting pileup summaries
process GetPileupSummaries {
    // Set maximum memory
    memory '40 GB'

    publishDir "${params.outdir}/summaries", mode: 'copy'

    input:
    path tumor_bam_sorted
    path normal_bam_sorted
    path exac
    val ID

    output:
    file("${tumor_bam_sorted.simpleName}.getpileupsummaries.table")
    file("${normal_bam_sorted.simpleName}.getpileupsummaries.table")
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
