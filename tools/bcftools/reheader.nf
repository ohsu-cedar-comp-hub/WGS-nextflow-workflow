#!/usr/bin/env nextflow

process REHEADER {
    
    publishDir "${params.outdir}/vcfs", mode: 'copy'
    container "${params.container_bcftools}"

    input:
    path filtered_vcf

    output: 
    file "${filtered_vcf.baseName}_reheader.vcf"
    
    script:

    """
    awk '/^##/ {print; next} {exit}' ${filtered_vcf} > temp_vcf_header.txt
    echo '##NextflowWGSPipelineVersion="v0.1.1 (2024-08-26) at https://github.com/ohsu-cedar-comp-hub/WGS-nextflow-workflow/releases/tag/v0.1.1\"' >> temp_vcf_header.txt
    awk '/^#CHROM/ {print; exit}' ${filtered_vcf} >> temp_vcf_header.txt
    bcftools reheader --header temp_vcf_header.txt ${filtered_vcf} > ${filtered_vcf.baseName}_reheader.vcf
    """
}