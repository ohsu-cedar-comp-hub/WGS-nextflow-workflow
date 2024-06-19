#!/usr/bin/env nextflow
params.vcf_files = "*.unfiltered.vcf"
params.outDir = "results"
Channel
    .fromFilePairs(params.vcf_files, size: 1, flat: true)
    .set { unfiltered_ch }

// concatenate VCF files
process CONCAT {
    publishDir path: "${params.outdir}/svc/concat", mode: 'copy'
    container "${params.container_bcftools}"
    
    tag { sample_id }
    
    input:
        tuple val(sample_id), path(vcf_files) from unfiltered_ch
    
    output:
        tuple val(sample_id), path("*.concat.vcf") into concatenated_ch
    
    script:
    """
    #!/bin/bash
    set -e
    bcftools concat -a ${vcf_files} -Ov -o ${sample_id}.concat.vcf
    """
}
// sort catted VCF files
process sort {
    publishDir params.outDir, mode: 'copy'
    tag { sample_id }
    input:
        tuple val(sample_id), path(concat_file) from concatenated_ch
    output:
        tuple val(sample_id), path("*.sorted.vcf") into sorted_ch
    script:
    """
    #!/bin/bash
    set -e
    bcftools sort ${concat_file} -Ov -o ${sample_id}.sorted.vcf
    """
}
// compress and index sorted vcf files
process compress {
    publishDir params.outDir, mode: 'copy'
    tag { sample_id }
    input:
        tuple val(sample_id), path(sorted_file) from sorted_ch
    output:
        tuple val(sample_id), path(".vcf.gz"), path(".vcf.gz.tbi") into compressed_ch
    script:
    """
    #!/bin/bash
    set -e
    # Compress
    bgzip -c ${sorted_file} > ${sample_id}.vcf.gz
    # Index
    bcftools index -t ${sample_id}.vcf.gz
    """
}