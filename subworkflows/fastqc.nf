// Define the process for running FastQC
process fastQC {
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    path Read1
    path Read2

    output:
    file("${Read1.baseName}_qc_report.html")
    file("${Read2.baseName}_qc_report.html")
    
    script:
    """
     fastqc ${Read1} > ${Read1.baseName}_qc_report.html 2>&1
     fastqc ${Read2} > ${Read2.baseName}_qc_report.html 2>&1
    """
}
