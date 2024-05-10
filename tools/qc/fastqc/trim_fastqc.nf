#!/usr/bin/env nextflow

// Define the process for running FastQC
process FastQC {
    publishDir "${params.outdir}/trimfastqc", mode: 'copy'

    input:
    path trim_read1
    path trim_read2
    val id
    
    script:
    """
     /usr/local/FastQC/fastqc -o ${params.outdir}/trimfastqc ${trim_read1}
     /usr/local/FastQC/fastqc -o ${params.outdir}/trimfastqc ${trim_read2} 
    """
}
