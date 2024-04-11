#!/usr/bin/env nextflow

workflow {
    // Run MultiQC on the specified output directory
    fastqc_read1=file(params.fastqc_read1)
    fastqc_read2=file(params.fastqc_read2)
    trim_fastqc_read1=file(params.trim_fastqc_read1)
    trim_fastqc_read2=file(params.trim_fastqc_read2)

    multiQC(fastqc_read1, fastqc_read2, trim_fastqc_read1, trim_fastqc_read2)
}
