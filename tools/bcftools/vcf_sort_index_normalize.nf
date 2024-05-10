#!/usr/bin/env nextflow << NOTE: This script can be used when input are not splitting by chromosome 

process vcfSortIndexNormalize {
    // Set output directory
    publishDir "${params.outdir}/svc/sort_index", mode: 'copy'

    // define input and output paramaters
    input: 
    path unfiltered_vcf
 

    output: 
    file("${unfiltered_vcf.baseName}_sorted.vcf.gz")
    file("${unfiltered_vcf.baseName}_normalize.vcf.gz") 

    script: 
    """
    // sort  command
    bcftools sort ${unfiltered_vcf} -Oz -o ${unfiltered_vcf.baseName}_sorted.vcf.gz
    // index command
    bcftools index -t ${unfiltered_vcf.baseName}_sorted.vcf.gz
    // normalize command
    bcftools norm -m -any ${unfiltered_vcf.baseName}_sorted.vcf.gz -o ${unfiltered_vcf.baseName}_normalize.vcf.gz
    """
}
