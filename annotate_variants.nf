#!/usr/bin/env nextflow

process Annotate_Variants {
    publishDir "${params.outdir}/svc/annotated_variants", mode: 'copy'

    input:
    path svc_vcf

    output: 
    file "${svc_vcf.baseName}_annotated_variants.vcf"

    script:
    """
    snpEff GRCh38.86 ${svc_vcf} -Xmx8g -cancer > ${svc_vcf.baseName}_annotated_variants.vcf
    """

}

workflow {
    // specify input sorted bam
    input_file=file(params.svc_vcf)

    //run MarkDuplications
    Annotate_Variants (input_file)

}
