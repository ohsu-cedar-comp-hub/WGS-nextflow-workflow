#!/usr/bin/env nextflow

process BGZIP {
    
    publishDir "${params.outdir}/svc", mode: 'copy'

    container "${params.container_bcftools}"

    input:
    path split_vcf

    output:
    path ("${split_vcf}.gz"), emit: vcf
    path ("${split_vcf}.gz.tbi"), emit: index // bcftools concat requires an index to accompany the bgzipped vcf

    script:
    """
    bgzip -c ${split_vcf} > ${split_vcf}.gz
    bcftools index -t ${split_vcf}.gz
    """
}

process PREPAREVCF {
    debug true
    
    publishDir "${params.outdir}/svc/sort_index", mode: 'copy'
    
    container "${params.container_bcftools}"

    input:
    path split_vcfs
    path split_vcfs_index
    val sample_id

    output:
    path("${sample_id}_unfiltered.vcf.gz"), emit: unfiltered
    path("${sample_id}_unfiltered_sorted.vcf.gz"), emit: sorted
    path("${sample_id}_unfiltered_sorted.vcf.gz.tbi"), emit: index
    path("${sample_id}_normalized.vcf.gz"), emit: normalized

    script:
    """
    bcftools concat -a ${split_vcfs.join(' ')} -o ${sample_id}_unfiltered.vcf.gz
    bcftools sort -Oz -o ${sample_id}_unfiltered_sorted.vcf.gz
    bcftools index -t ${sample_id}_unfiltered_sorted.vcf.gz
    bcftools norm -m -any ${sample_id}_unfiltered_sorted.vcf.gz -o ${sample_id}_normalized.vcf.gz
    """
}

