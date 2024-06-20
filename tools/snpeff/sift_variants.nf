#!/usr/bin/env nextflow

process SiftVariants {
    publishDir "${params.outdir}/svc/annotated_variants", mode: 'copy'
    // Set maximum memory
    memory '40 GB'

    input:
    path filtered_vcf
    val ID

    output: 
    file "${filtered_vcf.baseName}_passed.vcf"

    // GEN[0] is the germline
    // GEN[1] is the tumor
    // AD[0] is allelic depth of the reference
    // AD[1] is allelic depth of the alt
    
    script:
    """
    cat ${filtered_vcf} | java -Xmx8g -jar /usr/src/app/snpEff/SnpSift.jar filter "( ((GEN[0].AD[0] >= 3) | (GEN[0].AD[1] >= 4)) | ((GEN[1].AD[0] >= 3) | (GEN[1].AD[1] >= 4)) ) & ( ( na FILTER ) | (FILTER = 'PASS') )" > ${filtered_vcf.baseName}_passed.vcf
    """

}