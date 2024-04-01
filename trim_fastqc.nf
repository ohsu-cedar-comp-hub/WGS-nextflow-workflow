#!/usr/bin/env nextflow

// Define the process for running FastQC
process fastQC {
    publishDir "${params.outdir}/trimfastqc", mode: 'copy'

    input:
    path trim_read1
    path trim_read2
    
    script:
    """
     fastqc ${trim_read1} -o ${params.outdir}/trimfastqc
     fastqc ${trim_read2} -o ${params.outdir}/trimfastqc 
    """
}

workflow {
    trim_read1=file(params.trim_read1)
    trim_read2=file(params.trim_read2)

    quality_check_results = fastQC(trim_read1, trim_read2)
}
