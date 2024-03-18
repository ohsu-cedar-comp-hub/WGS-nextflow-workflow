#!/usr/bin/env nextflow

// Define the process for running FastQC
process fastQC {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file Read1 
    file Read2

    output:
    file("*_R1_fastqc.html")
    file("*_R2_fastqc.html")
    
    script:
    """
    fastqc ${Read1} -o ${params.outdir}
    fastqc ${Read2} -o ${params.outdir}
    """
}

workflow {
    // Run FastQC for each specified fastq file
    fastQC(params.Read1, params.Read2)
}