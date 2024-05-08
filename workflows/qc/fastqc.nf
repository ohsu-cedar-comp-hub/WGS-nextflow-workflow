#!/usr/bin/env nextflow

include {fastQC} from '../../tools/qc/fastqc/fastqc.nf'

// Define the workflow for fastqc

workflow {

    fastQC(file(params.read1), file(params.read2), "test")
}