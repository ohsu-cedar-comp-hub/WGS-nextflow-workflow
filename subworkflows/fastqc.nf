// Define the process for running FastQC
process fastQC {
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    path read1
    path read2

    output:
    file("${read1.baseName}_qc_report.html")
    file("${read2.baseName}_qc_report.html")
    
    script:
    """
     fastqc ${read1} > ${read1.baseName}_qc_report.html 2>&1
     fastqc ${read2} > ${read2.baseName}_qc_report.html 2>&1
    """
}