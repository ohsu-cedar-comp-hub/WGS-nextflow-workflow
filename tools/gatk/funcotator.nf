#!/usr/bin/env nextflow

process FUNCOTATOR {

    container "${params.container_gatk}"

    publishDir "${params.outdir}/vcfs", mode: 'copy'

    input:
    path vcf
    path mutect_idx
    path mutect_idx_fai
    path mutect_idx_dict
    path funcotator_data
    val sample_id

    output:
    path("${vcf.basename}_FUNCOTATED.vcf")

    script:
    """
    gatk Funcotator \
    --variant ${vcf} \
    --reference ${mutect_idx} \
    --ref-version hg38 \
    --data-sources-path ${funcotator_data} \
    --output ${vcf.basename}_FUNCOTATED.vcf \
    --output-file-format VCF
    """
}
