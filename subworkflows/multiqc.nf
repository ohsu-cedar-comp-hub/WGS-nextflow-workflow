process multiQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    params.outdir
    output:
    file("multiqc_report.html")

    script:
    """
    multiqc ${params.outdir} -o ${params.outdir}/multiqc
    """
}
