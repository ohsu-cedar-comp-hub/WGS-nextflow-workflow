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

    output:
    file("${tumor_bam_sorted.baseName}.getpileupsummaries.table")
    file("${normal_bam_sorted.baseName}.getpileupsummaries.table")
    script:
    """
    //script for tumor_bam
    gatk GetPileupSummaries \\
    -I ${tumor_bam_sorted} \\
    -V ${exac} \\
    -L ${exac} \\
    -O ${tumor_bam_sorted.baseName}.getpileupsummaries.table

    //script for normal_bam
    gatk GetPileupSummaries \\
    -I ${normal_bam_sorted} \\
    -V ${exac} \\
    -L ${exac} \\
    -O ${normal_bam_sorted.baseName}.getpileupsummaries.table
    """
}
