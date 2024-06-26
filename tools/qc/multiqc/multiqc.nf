#!/usr/bin/env nextflow

//Define process for multiQC
process MultiQC {
    // Set output directory for multiQC results
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    // Define input paramaters from fastqc runs
    input:
    path fastqc_read1
    path fastqc_read2
    path trim_fastqc_read1
    path trim_fastqc_read2
    val id

    output:
    file("${fastqc_read1.baseName}_multiqc_report.html")
    file("${fastqc_read2.baseName}_multiqc_report.html")
    file("${trim_fastqc_read1.baseName}_multiqc_report.html")
    file("${trim_fastqc_read2.baseName}_multiqc_report.html")

    //multiqc command
    script:
    """
    multiqc ${fastqc_read1} -o ${fastqc_read1.baseName}_multiqc_report.html
    multiqc ${fastqc_read2} -o ${fastqc_read2.baseName}_multiqc_report.html
    multiqc ${trim_fastqc_read1} -o ${trim_fastqc_read1.baseName}_multiqc_report.html
    multiqc ${trim_fastqc_read2} -o ${trim_fastqc_read2.baseName}_multiqc_report.html
    """
}
