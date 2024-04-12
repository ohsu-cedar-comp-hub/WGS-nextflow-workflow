#!/usr/bin/env nextflow

include { bwaMem2Alignment } from '../tools/bwa/paired_align.nf'

// Define the workflow
workflow {
    // Run BWA-MEM2 alignment for each read file
    bwaMem2Alignment(file(params.read1), file(params.read2), file(params.idx), "test")
}
