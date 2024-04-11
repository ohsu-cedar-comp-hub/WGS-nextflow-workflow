#!/usr/bin/env nextflow

// define process GATK4 CalculateContamination
process CalculateContamination {
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
        -I ${params.tumor_pileups_table} \\
        --matched ${params.normal_pileups_table} \\
        -O ${tumor_pileups_table.baseName}_contamination_table \\
        -tumor-segmentation ${tumor_pileups_table.baseName}_segmentation_table
    """
  // Define the workflow
workflow {
    // define input paramaters
    input_tumor_table=file(params.tumor_pileups_table)
    input_matched_table=file(params.normal_pipeups_table)
    // Run the workflow
    CalculateContamination(input_tumor_table, input_matched_table)
}
