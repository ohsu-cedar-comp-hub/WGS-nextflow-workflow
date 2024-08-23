#!/usr/bin/env nextflow

process SNPEFF {

    publishDir "${params.outdir}/vcfs", mode: 'copy', pattern: '*.vcf'
    container "${params.container_snpeff}"

    input:
    path filtered_vcf
    val sample_id

    output: 
    path("${filtered_vcf.baseName}_SNPANNOTATED.vcf"), emit: vcf

    script:
    """
    java -Xmx8g -jar /usr/src/app/snpEff/snpEff.jar GRCh38.86 ${filtered_vcf} -cancer > ${filtered_vcf.baseName}_SNPANNOTATED.vcf
    """

}