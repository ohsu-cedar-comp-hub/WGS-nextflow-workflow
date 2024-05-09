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
     /usr/local/FastQC/fastqc -o ${params.outdir}/fastqc ${read1}  
     /usr/local/FastQC/fastqc -o ${params.outdir}/fastqc ${read2} 
    """
}
