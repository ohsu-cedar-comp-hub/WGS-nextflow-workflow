#!/usr/bin/env nextflow

include {FastQC} from '../../tools/qc/fastqc/trim_fastqc.nf'	

// Define the workflow for fastqc

workflow {
    FastQC(file(params.trim_read1), file(params.trim_read2), "test")
}
