#!/usr/bin/env nextflow

process SNPSIFT {
    
    publishDir "${s3outdir}/annotated", mode: 'copy'
    container "${params.container_snpeff}"

    input:
    path filtered_vcf
    val sample_id

    output: 
    file "${sample_id}_annotated_PASSED.vcf"

    // GEN[0] is the germline
    // GEN[1] is the tumor
    // AD[0] is allelic depth of the reference
    // AD[1] is allelic depth of the alt
    
    script:
    """
    cat ${filtered_vcf} | java -Xmx8g -jar /usr/src/app/snpEff/SnpSift.jar filter "( ((GEN[0].AD[0] >= 3) | (GEN[0].AD[1] >= 4)) | ((GEN[1].AD[0] >= 3) | (GEN[1].AD[1] >= 4)) ) & ( ( na FILTER ) | (FILTER = 'PASS') )" > ${sample_id}_annotated_PASSED.vcf
    """

}