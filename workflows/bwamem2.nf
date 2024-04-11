#!/usr/bin/env nextflow

// Define the workflow
workflow {
    // Define input parameters
    trim_read1=file(params.trim_read1)
    trim_read2=file(params.trim_read2)
    idx=file(params.idx)

    // Run BWA-MEM2 alignment for each read file
    bwaMem2Alignment(trim_read1, trim_read2, idx)
}
