#!/usr/bin/env nextflow

include { BwaMem2Alignment } from '../../tools/bwa/bwamem2.nf'

// Define the workflow
workflow {
    // Run BWA-MEM2 alignment for each read file
    BwaMem2Alignment(file(params.trim_read1), file(params.trim_read2), file(params.idx), params.ID)
}
