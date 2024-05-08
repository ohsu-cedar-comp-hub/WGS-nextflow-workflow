#!/usr/bin/env nextflow

// Define the process for running FastQC
process fastQC {
    publishDir "${params.outdir}/trimfastqc", mode: 'copy'

    input:
    path trim_read1
    path trim_read2
    val id
    
    script:
    """
     fastqc -o ${params.outdir}/trimfastqc ${trim_read1}
     fastqc -o ${params.outdir}/trimfastqc ${trim_read2} 
    """
}
