#!/usr/bin/env nextflow

process vcfSortIndexNormalize {
    // Set output directory
    publishDir "${params.outdir}/svc/sort_index", mode: 'copy'

    // define input and output paramaters
    input: 
    path unfiltered_vcf
    val ID
 

    output: 
    file("${unfiltered_vcf.baseName}_sorted.vcf.gz")
    file("${unfiltered_vcf.baseName}_normalize.vcf.gz") 

    script: 
    """
    bcftools sort ${unfiltered_vcf} -Oz -o ${unfiltered_vcf.baseName}_sorted.vcf.gz
    bcftools index -t ${unfiltered_vcf.baseName}_sorted.vcf.gz
    bcftools norm -m -any ${unfiltered_vcf.baseName}_sorted.vcf.gz -o ${unfiltered_vcf.baseName}_normalize.vcf.gz
    """
}
