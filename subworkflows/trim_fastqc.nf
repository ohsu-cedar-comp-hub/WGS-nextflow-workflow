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
