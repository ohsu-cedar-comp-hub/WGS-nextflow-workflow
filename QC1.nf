#!/usr/bin/env nextflow


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

workflow {
    input_file=file(params.infile)

    quality_check_results = qualityControl(input_file)
}

