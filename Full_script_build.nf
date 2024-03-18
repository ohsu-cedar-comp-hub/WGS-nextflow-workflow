#!/usr/bin/env nextflow

process simpleMD5 {
    publishDir "$params.outdir"
    input:
        path 'infile'

    output:
        path 'out_file'

    script:
        """
        md5sum $infile > out_file
        """
}

process simpleWC {
    publishDir "$params.outdir"
    input:
        path sourcefile

    output:
        path 'out_wc'

    script:
        """
        wc $sourcefile > out_wc
        """
}

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
    a = simpleMD5(file(params.infile))
    simpleWC(a)

    input_file=file(params.infile)
    quality_check_results = qualityControl(input_file)
}



