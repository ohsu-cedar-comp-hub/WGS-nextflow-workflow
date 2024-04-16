#!/usr/bin/env nextflow

// define process GATK4 CalculateContamination
process CalculateContamination {
    // Set maximum memory
    memory '40 GB'

    // Set output directory
    publishDir "${params.outdir}/tables", mode: 'copy'
  
  // define input and output
    input:
    path tumor_pileups_table
    path normal_pileups_table

    output:
    path("${tumor_pileups_table.baseName"}_contamination_table)
    path("${tumor_pileups_table.baseName"}_segmentation_table)

  // calculatContamination command
    script:
    """
        gatk4 CalculateContamination \\
        -I ${tumor_pileups_table} \\
        --matched ${normal_pileups_table} \\
        -O ${tumor_pileups_table.baseName}_contamination_table \\
        -tumor-segmentation ${tumor_pileups_table.baseName}_segmentation_table
