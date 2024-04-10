#!/usr/bin/env nextflow

// Define the workflow
workflow {
    // define input paramaters for Mutect2
    input_file=file(params.bam_sorted)
    idx=file(params.mutect_idx)
    pon=file(params.pon)
    // Run Mutect2
    mutect2(input_file, idx, pon)
}
