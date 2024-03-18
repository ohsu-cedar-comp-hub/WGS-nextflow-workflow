#!/usr/bin/env nextflow

// Define the process for running FastQC
process qualityControl {
    publishDir "${params.outdir}"
    input:
    path infile

    output:
    file("${infile.baseName}_qc_report.html")

    script:
    """
    fastqc ${infile} > ${infile.baseName}_qc_report.html 2>&1
    """
}

// Define the process for running MultiQC
process multiQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    file('*_qc_report.html') 

    output:
    file("multiqc_report.html")

    script:
    """
    multiqc ${params.outdir}/fastqc -o ${params.outdir}/multiqc
    """
}

workflow {
    // Run FastQC for each specified fastq file
    input_file=file(params.infile)
    quality_check_results = qualityControl(input_file)

    // Run MultiQC on the FastQC output
    multiQC(quality_check_results)
}
