#!/usr/bin/env nextflow

process PREPAREVCF {
    publishDir "${params.outdir}/svc/sort_index", mode: 'copy'
    
    container "${params.container_bcftools}"

    input:
    path split_vcfs
    val sample_id

    output:
    path("${sample_id}_unfiltered_all.vcf.bgz")
    path("${sample_id}_unfiltered_all.vcf.bgz.tbi")

    script:
    """
    files_to_concat=""
    for chr in {1..22} X; do
        files_to_concat+=" ${unfiltered_vcf.baseName}_chr\${chr}_unfiltered.vcf.bgz"
    done

    bcftools concat -a $split_vcfs -o ${sample_id}_unfiltered_all.vcf.bgz
    bcftools index -t ${sample_id}_unfiltered_all.vcf.bgz
    """
}