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
    bcftools filter -i "(FORMAT/DP[0] > 6) && (FORMAT/DP[1] > 10) && (FORMAT/AD[0:0] > 5) && (FORMAT/AD[0:1] > 3) && (FORMAT/AD[1:0] > 5) && (FORMAT/AF[0:1] > 0.05) && (INFO/MMQ[0] > 40) && (INFO/MMQ[1] > 40)" -o ${filtered_vcf.baseName}_bcffilter.vcf ${filtered_vcf}
    """
}