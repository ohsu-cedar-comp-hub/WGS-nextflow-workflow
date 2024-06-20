#!/usr/bin/env nextflow

process BGZIP {

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
    
    publishDir path: "${params.outdir}/svc/sort_index", mode: 'copy', pattern: "*_unfiltered_normalized_sorted.vcf.gz"
    
    container "${params.container_bcftools}"

    input:
    path split_vcfs
    path split_vcfs_index
    val sample_id
    path mutect_idx
    path mutect_idx_fai
    path mutect_idx_dict
    
    output:
    path("${sample_id}_unfiltered.vcf.gz"), emit: vcf
    path("${sample_id}_unfiltered_normalized_sorted.vcf.gz"), emit: normalized
    path("${sample_id}_unfiltered_normalized_sorted.vcf.gz.tbi"), emit: index

    script:
    """
    bcftools concat -a ${split_vcfs.join(' ')} -o ${sample_id}_unfiltered.vcf.gz
    bcftools sort -Oz ${sample_id}_unfiltered.vcf.gz -o ${sample_id}_unfiltered_sorted.vcf.gz 
    bcftools index -t ${sample_id}_unfiltered_sorted.vcf.gz
    bcftools norm -Oz -f ${mutect_idx} ${sample_id}_unfiltered_sorted.vcf.gz -o ${sample_id}_unfiltered_normalized.vcf.gz
    bcftools sort -Oz ${sample_id}_unfiltered_normalized.vcf.gz -o ${sample_id}_unfiltered_normalized_sorted.vcf.gz
    bcftools index -t ${sample_id}_unfiltered_normalized_sorted.vcf.gz
    """
}



