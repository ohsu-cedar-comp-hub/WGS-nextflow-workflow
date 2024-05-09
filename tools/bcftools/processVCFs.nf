#!/usr/bin/env nextflow

// Define the process to prepare VCFs
process prepareVCFs {
    // Set output directory
    publishDir "${params.outdir}/svc/sort_index", mode: 'copy'
    
    input:
    path unfiltered_vcf

    output:
    path("${unfiltered_vcf.baseName}_unfiltered_all.vcf.bgz")
    path("${unfiltered_vcf.baseName}_unfiltered_all.vcf.bgz.tbi")

    script:
    """
    files_to_concat=""
    for chr in {1..22} X; do
        files_to_concat+=" ${unfiltered_vcf.baseName}_chr\${chr}_unfiltered.vcf.bgz"
    done

    bcftools concat -a \$files_to_concat -o ${unfiltered_vcf.baseName}_unfiltered_all.vcf.bgz
    bcftools index -t ${unfiltered_vcf.baseName}_unfiltered_all.vcf.bgz
    """
}