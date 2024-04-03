#!/usr/bin/env nextflow

process Annotate_Variants {
    publishDir "${params.outdir}/svc/annotated_variants", mode: 'copy'

    input:
    path svc_vcf

    output: 
    file "${svc_vcf.baseName}_annotated_variants.vcf"

    script:
    """
    snpEff GRCh38.86 ${svc_vcf} -Xmx4g > ${svc_vcf.baseName}_annotated_variants.vcf
    """

}
