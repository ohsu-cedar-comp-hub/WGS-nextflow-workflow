#!/usr/bin/env nextflow

// Process for getting pileup summaries
process GETPILEUPSUMMARIES {
    maxForks 3
    container "${params.container_gatk}"

    publishDir "${params.outdir}/intermediates", mode: 'copy'

    input:
    path bam_file
    path exac

    output:
    path ("*${params.tumor}*.getpileupsummaries.table"), emit: tumor, optional: true
    path ("*${params.normal}*.getpileupsummaries.table"), emit: normal, optional: true
    
    script:
    """
    gatk GetPileupSummaries \\
    -I ${bam_file} \\
    -V ${params.exac} \\
    -L ${params.exac} \\
    -O ${bam_file.baseName}.getpileupsummaries.table
    """
}
