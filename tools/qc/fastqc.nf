#!/usr/bin/env nextflow

// Define the process for running FastQC
process fastQC {
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    path read1
    path read2
    val ID
    
    script:
    """
     fastqc -o ${params.outdir}/fastqc ${read1}  
     fastqc -o ${params.outdir}/fastqc ${read2} 
    """
}
