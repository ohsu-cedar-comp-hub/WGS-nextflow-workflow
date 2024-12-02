#!/usr/bin/env nextflow

// Process for getting pileup summaries
process GETPILEUPSUMMARIES {
    maxForks 3
    container "${params.container_gatk}"

    publishDir "${params.outdir}/intermediates", mode: 'copy'

    input:
    path bam_file
    path bai_file
    path exac

    output:
    path ("*${params.tumor_namepattern}*.getpileupsummaries.table"), emit: tumor, optional: true
    path ("*${params.normal_namepattern}*.getpileupsummaries.table"), emit: normal, optional: true
    
    script:
    """
    gatk GetPileupSummaries \\
    -I ${bam_file} \\
    -V ${params.exac} \\
    -L ${params.exac} \\
    -O ${bam_file.baseName}.getpileupsummaries.table
    """
}

process TUMORONLYGETPILEUPSUMMARIES {
    maxForks 2
    container "${params.container_gatk}"

    publishDir "${params.outdir}/intermediates", mode: 'copy'

    input:
    path bam_file
    path bai_file
    path exac

    output:
    path ("*.getpileupsummaries.table")
    
    script:
    """
    gatk GetPileupSummaries \\
    -I ${bam_file} \\
    -V ${params.exac} \\
    -L ${params.exac} \\
    -O ${bam_file.baseName}.getpileupsummaries.table
    """
}