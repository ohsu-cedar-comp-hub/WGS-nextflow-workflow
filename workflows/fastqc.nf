#!/usr/bin/env nextflow

workflow {
    read1=file(params.read1)
    read2=file(params.read2)

    quality_check_results = fastQC(read1, read2)
}
