#!/usr/bin/env nextflow

include { bwaMem2Alignment } from '../tools/bwa/bwamem2.nf'

// Define the workflow
workflow {
    // Run BWA-MEM2 alignment for each read file
    bwaMem2Alignment(file(params.trim_read1), file(params.trim_read2), file(params.idx), "test")
}
