#!/usr/bin/env nextflow

workflow {
    trim_read1=file(params.trim_read1)
    trim_read2=file(params.trim_read2)

    quality_check_results = fastQC(trim_read1, trim_read2)
}
