#!/usr/bin/env nextflow

process SiftVariants {
    publishDir "${params.outdir}/svc/annotated_variants", mode: 'copy'
    // Set maximum memory
    memory '40 GB'

    input:
    path filtered_vcf
    val ID

    output: 
    file "${filtered_vcf.baseName}_passed.vcf"

    script:
    """
    cat ${filtered_vcf} | java -Xmx8g -jar /usr/src/app/snpEff/SnpSift.jar filter "( na FILTER ) | (FILTER = 'PASS')" > ${filtered_vcf.baseName}_passed.vcf
    """

}