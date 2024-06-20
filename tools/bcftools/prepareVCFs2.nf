#!/usr/bin/env nextflow

// compress and index sorted vcf files
process COMPRESS {
    
    publishDir path: "${params.outdir}/svc/compress", mode: 'copy'
    container "${params.container_bcftools}"
    
    input:
    path vcf
    val sample_id
    
    output:
    path("${vcf}.gz"), emit: vcf
    path("${vcf}.gz.tbi"), emit: tbi
    
    script:
    """
    #!/bin/bash
    set -e
    # Compress
    bgzip -c ${vcf} > ${vcf}.gz
    # Index
    bcftools index -t ${vcf}.gz
    """
}

// concatenate VCF files
process CONCAT {
    publishDir path: "${params.outdir}/svc/concat", mode: 'copy'
    container "${params.container_bcftools}"
    
    input:
    path split_vcfs
    path split_vcfs_index
    val sample_id
    
    output:
    path("${sample_id}_concat.vcf")
    
    script:
    """
    #!/bin/bash
    set -e
    bcftools concat -a ${split_vcfs.join(' ')} -Ov -o ${sample_id}_concat.vcf
    """
}

// sort catted VCF files
process SORT {
    
    publishDir path: "${params.outdir}/svc/sort", mode: 'copy'
    container "${params.container_bcftools}"

    input:
    path concat_vcf
    val sample_id
    
    output:
    path("${sample_id}_sorted.vcf")
    
    script:
    """
    #!/bin/bash
    set -e
    bcftools sort ${concat_vcf} -Ov -o ${sample_id}_sorted.vcf
    """
}

process NORM {

    publishDir path: "${params.outdir}/svc/normalize", mode: 'copy'
    container "${params.container_bcftools}"

    input:
    path compress_vcf
    path compress_vcf_tbi
    val sample_id

    output: 
    path("${sample_id}_normalized.vcf.gz"), emit: normalized
    path("${sample_id}_normalized.vcf.gz.tbi"), emit: index
    
    script:
    """
    #!/bin/bash
    set -e
    bcftools norm -Oz -m -any ${compress_vcf} -o ${sample_id}_normalized.vcf.gz
    bcftools index -t ${sample_id}_normalized.vcf.gz
    """
}