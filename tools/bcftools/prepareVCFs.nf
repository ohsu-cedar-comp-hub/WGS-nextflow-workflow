#!/usr/bin/env nextflow

process BGZIP {
    
    container "${params.container_samtools}"

    input:
    path split_vcf

    output:
    path ("${split_vcf}.gz"), emit: vcf

    script:
    """
    bgzip ${split_vcf}
    """
}

process PREPAREVCF {
    debug true
    
    publishDir "${params.outdir}/svc/sort_index", mode: 'copy'
    
    container "${params.container_bcftools}"

    input:
    val split_vcfs
    val sample_id

    output:
    path("${sample_id}_unfiltered.vcf.gz"), emit: unfiltered
    path("${sample_id}_unfiltered_sorted.vcf.gz"), emit: sorted
    path("${sample_id}_unfiltered_sorted.vcf.gz.tbi"), emit: index
    path("${sample_id}_normalized.vcf.gz"), emit: normalized

    script:
    """
    echo ${split_vcfs.join(' ')}
    bcftools concat -a ${split_vcfs.join(' ')} -o ${sample_id}_unfiltered.vcf.gz
    bcftools sort -Oz -o ${sample_id}_unfiltered_sorted.vcf.gz
    bcftools index -t ${sample_id}_unfiltered_sorted.vcf.gz
    bcftools norm -m -any ${sample_id}_unfiltered_sorted.vcf.gz -o ${sample_id}_normalized.vcf.gz
    """
}

