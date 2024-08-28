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
    echo '##NextflowWGSPipelineVersion="${params.release} (${params.releasedate}) at ${params.githublink}\"' >> temp_vcf_header.txt
    awk '/^#CHROM/ {print; exit}' ${filtered_vcf} >> temp_vcf_header.txt
    bcftools reheader --header temp_vcf_header.txt ${filtered_vcf} > ${filtered_vcf.baseName}_reheader.vcf
    """
}