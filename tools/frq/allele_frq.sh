#!/usr/bin/bash

# Example usage:
# bash allele_frq.sh src/G_H_507472.vcf vcftools/G_H_507472
# bash allele_frq.sh src/T_H_50472.vcf vcftools/T_H_50472

VCF=${1}
OUT=${2}

# Allele freq. exclude sites with more than 2 alleles (creates .frq)
echo 'Running allele freq'
vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2

# site quality (creates .lqual)
echo 'Running site qual'
vcftools --gzvcf $VCF --site-quality --out $OUT

# mean depth per site (creates *.ldepth.mean)
echo 'Running mean depth per site'
vcftools --gzvcf $SUBSET_VCF --site-mean-depth --out $OUT
