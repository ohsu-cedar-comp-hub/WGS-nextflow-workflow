#!/usr/bin/env nextflow

// Define the process for running FastQC
process fastQC {
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    path read1
    path read2
    
    script:
    """
     fastqc -o ${params.outdir}/fastqc ${read1}  
     fastqc -o ${params.outdir}/fastqc ${read2} 
    """
}

workflow {
    read1=file(params.read1)
    read2=file(params.read2)

    quality_check_results = fastQC(read1, read2)
}
