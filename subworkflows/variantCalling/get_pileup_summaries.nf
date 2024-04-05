#!/usr/bin/env nextflow

// Process for getting pileup summaries
process GetPileupSummaries {
    publishDir "${params.outdir}/summaries", mode: 'copy'

    input:
    path bam_sorted
    path exac

    output:
    file("${bam_sorted.baseName}.getpileupsummaries.table")

    script:
    """
    gatk GetPileupSummaries \\
    -I ${params.bam_sorted} \\
    -V ${params.exac} \\
    -L ${params.exac} \\
    -O ${bam_sorted.baseName}.getpileupsummaries.table
    """
}
