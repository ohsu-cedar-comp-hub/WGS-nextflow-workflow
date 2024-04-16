#!/usr/bin/env nextflow

process Annotate_Variants {
    publishDir "${params.outdir}/svc/annotated_variants", mode: 'copy'
    // Set maximum memory
    memory '40 GB'

    input:
    path filtered_vcf

    output: 
    file "${filtered_vcf.baseName}_annotated_variants.vcf"

    script:
    """
    snpEff GRCh38.86 ${filtered_vcf} -Xmx8g -cancer > ${filtered_vcf.baseName}_annotated_variants.vcf
    """

}
