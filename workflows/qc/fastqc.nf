#!/usr/bin/env nextflow

include {FastQC} from '../../tools/qc/fastqc/fastqc.nf'

// Define the workflow for fastqc

workflow {

    FastQC(file(params.read1), file(params.read2), params.ID)
}
