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
    bcftools filter -i "(INFO/GERMQ < 20 || INFO/NLOD < 2.2) && (FORMAT/DP[0] > 10) && (FORMAT/DP[1] > 6) && (FORMAT/AD[0:0] > 5) && (FORMAT/AD[1:0] > 5) && (FORMAT/AD[1:1] > 5) && (FORMAT/AF[1:0] > 0.05) && (INFO/MMQ[0] > 40) && (INFO/MMQ[1] > 40) && (INFO/gnomad_AF < 0.01) && (INFO/gnomad_AF_afr < 0.01)" -o ${filtered_vcf.baseName}_bcffilter.vcf ${filtered_vcf}
    """
}