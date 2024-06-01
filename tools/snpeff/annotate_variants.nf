#!/usr/bin/env nextflow

process AnnotateVariants {
    
    container "/home/groups/CEDAR/lancasru/WGS_COH_NF/config_sif/snpeff_conda.sif"
    conda "/home/groups/CEDAR/lancasru/anaconda3/envs/nextflow_env"

    publishDir "${params.outdir}/svc/annotated_variants", mode: 'copy'
    // Set maximum memory
    memory '40 GB'

    input:
    path filtered_vcf
    val ID

    output: 
    file "${filtered_vcf.baseName}_annotated_variants.vcf"

    script:
    """
    java -Xmx8g -jar /usr/src/app/snpEff/snpEff.jar GRCh38.86 ${filtered_vcf} -Xmx8g -cancer > ${filtered_vcf.baseName}_annotated_variants.vcf
    """

}