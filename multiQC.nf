#!/usr/bin/env nextflow

// Define the process for running FastQC
process fastQC {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file Read1
    file Read2

    output:
    file("${Read1.baseName}_R1_fastqc.html")
    file("${Read2.baseName}_R2_fastqc.html")
    
    script:
    """
    fastqc ${Read1} -o ${params.outdir}
    fastqc ${Read2} -o ${params.outdir}
    """
}

// Define the process for running MultiQC
process multiQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    file("${Read1.baseName}_R1_fastqc.html")
    file("${Read2.baseName}_R2_fastqc.html")

    output:
    file("multiqc_report.html")

    script:
    """
    multiqc ${params.outdir}/fastqc -o ${params.outdir}/multiqc
    """
}

workflow {
// Run FastQC for each specified fastq file
    Read1=file(params.Read1)
    Read2=file(params.Read2)

    fastQC(Read1, Read2)

    // Run MultiQC on the FastQC output
    multiQC(quality_check_results)
}
