#!/usr/bin/env nextflow

process multiQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path fastqc_read1
    path fastqc_read2
    output:
    file("${fastqc_read1.baseName}_multiqc_report.html")

    script:
    """
    multiqc ${fastqc_read1} ${fastqc_read2} -o ${fastqc_read1.baseName}_multiqc_report.html
    """
}

