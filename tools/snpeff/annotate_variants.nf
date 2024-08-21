#!/usr/bin/env nextflow

process SNPEFF {

    publishDir "${params.outdir}/vcfs", mode: 'copy'
    container "${params.container_snpeff}"

    input:
    path filtered_vcf
    val sample_id

    output: 
    path("${sample_id}_PASSED_annotated.vcf")

    script:
    """
    java -Xmx8g -jar /usr/src/app/snpEff/snpEff.jar GRCh38.86 ${filtered_vcf} -cancer > ${sample_id}_PASSED_annotated.vcf
    """

}