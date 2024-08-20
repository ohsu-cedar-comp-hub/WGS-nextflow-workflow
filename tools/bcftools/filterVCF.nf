#!/usr/bin/env nextflow

process ADDFILTER {
    
    publishDir "${params.outdir}/svc", mode: 'copy'
    container "${params.container_bcftools}"

    input:
    path filtered_vcf
    val sample_id

    output: 
    file "${filtered_vcf.basename}_bcffilter.vcf"
    
    script:
    """
    bcftools filter -i "(FORMAT/DP[0] > 10) 
                && (FORMAT/DP[1] > 6) 
                && (FORMAT/AD[0:0] > 5) 
                && (FORMAT/AD[1:0] > 5) 
                && (FORMAT/AD[1:1] > 5)
                && (FORMAT/AF[1:0] > 0.05)
                && (INFO/MMQ[0] > 40)
                && (INFO/MMQ[1] > 40) " \
            -o ${filtered_vcf.basename}_bcffilter.vcf \
            ${filtered_vcf}
    """
}