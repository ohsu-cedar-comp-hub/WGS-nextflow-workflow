#!/usr/bin/env nextflow

process multiQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    params.outdir/fastqc
    output:
    file("multiqc_report.html")

    script:
    """
    multiqc ${params.outdir}/fastqc -o ${params.outdir}/multiqc
    """
}

workflow {
    // Run MultiQC on the specified output directory
    multiQC()

}
