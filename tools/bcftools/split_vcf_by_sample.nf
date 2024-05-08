#!/usr/bin/env nextflow

process SplitVcfBySample {

publishDir "${params.outdir}/svc/annotated_variants", mode: 'copy'

input:
path annotated_vcf

output:
file "${annotated_vcf.baseName}_normal.vcf"
file "${annotated_vcf.baseName}_tumor.vcf"

script: 
"""
# Filter sample names with a regular expression
normal_sample=\$(bcftools query -l ${annotated_vcf} | grep -E '^G')
tumor_sample=\$(bcftools query -l ${annotated_vcf} | grep -E '^T')

bcftools view -O v -s \$normal_sample ${annotated_vcf} -o ${annotated_vcf.baseName}_normal.vcf
bcftools view -O v -s \$tumor_sample ${annotated_vcf} -o ${annotated_vcf.baseName}_tumor.vcf

"""
}
