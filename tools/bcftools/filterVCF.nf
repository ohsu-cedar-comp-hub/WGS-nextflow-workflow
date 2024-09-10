#!/usr/bin/env nextflow

process ADDFILTER {
    
    publishDir "${params.outdir}/vcfs", mode: 'copy'
    container "${params.container_bcftools}"

    input:
    path filtered_vcf
    val sample_id

    output: 
    file "${filtered_vcf.baseName}_bcffilter.vcf"
    
    script:
    // the command has to be in one long line otherwise it isn't recorded in its entirety in the VCF header.
    """
    bcftools filter -i "${params.bcftools_filter}" -o ${filtered_vcf.baseName}_bcffilter.vcf ${filtered_vcf}
    """
}