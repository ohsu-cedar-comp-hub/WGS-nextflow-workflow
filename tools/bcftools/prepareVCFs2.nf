#!/usr/bin/env nextflow

// concatenate VCF files
process CONCAT {
    publishDir path: "${params.outdir}/svc/concat", mode: 'copy'
    container "${params.container_bcftools}"
    
    input:
    path split_vcfs
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

// compress and index sorted vcf files
process COMPRESS {
    
    publishDir path: "${params.outdir}/svc/compress", mode: 'copy'
    container "${params.container_bcftools}"
    
    input:
    path sorted_vcf
    val sample_id
    
    output:
    path("${sample_id}.vcf.gz"), emit: vcf
    path("${sample_id}.vcf.gz.tbi"), emit: tbi
    
    script:
    """
    #!/bin/bash
    set -e
    # Compress
    bgzip -c ${sorted_vcf} > ${sample_id}.vcf.gz
    # Index
    bcftools index -t ${sample_id}.vcf.gz
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