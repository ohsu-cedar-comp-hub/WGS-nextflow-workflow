#!/usr/bin/env nextflow

// define process GATK4 CalculateContamination
process CALCULATECONTAMINATION {

    container "${params.container_gatk}"

    publishDir "${params.outdir}/tables", mode: 'copy'
  
    input:
    path tumor_pileups_table
    path normal_pileups_table

    output:
    path ("${tumor_pileups_table.simpleName}_contamination_table"), emit: contamination
    path ("${tumor_pileups_table.simpleName}_segmentation_table"), emit: segment

    script:
    """
        gatk CalculateContamination \\
        -I ${tumor_pileups_table} \\
        --matched ${normal_pileups_table} \\
        -O ${tumor_pileups_table.simpleName}_contamination_table \\
        -tumor-segmentation ${tumor_pileups_table.simpleName}_segmentation_table
    """
}