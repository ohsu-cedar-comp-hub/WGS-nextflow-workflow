#!/usr/bin/env nextflow

process PASS {
    
    publishDir "${params.outdir}/vcfs", mode: 'copy'
    container "${params.container_snpeff}"

    input:
    path filtered_vcf
    val sample_id

    output: 
    file "${sample_id}_PASSED.vcf"
    
    script:
    """
    cat ${filtered_vcf} | java -Xmx8g -jar /usr/src/app/snpEff/SnpSift.jar filter "( ( na FILTER ) | (FILTER = 'PASS') )" > ${sample_id}_PASSED.vcf
    """

}