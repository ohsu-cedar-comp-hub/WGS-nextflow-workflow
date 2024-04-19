#!/usr/bin/bash

# Example usage:
# bash allele_frq.sh src/G_H_507472.vcf vcftools/G_H_507472
# bash allele_frq.sh src/T_H_50472.vcf vcftools/T_H_50472

VCF=${1}
OUT=${2}

# Allele freq. exclude sites with more than 2 alleles
echo 'Running allele freq'
vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2
