#!/usr/bin/env nextflow

// Define the workflow
workflow {
    // Define input parameters
    tumor_bam = file(params.bam_sorted)
    exac = file(params.exac)

    // Run the GetPileupSummaries process
    GetPileupSummaries(tumor_bam,exac)
}
