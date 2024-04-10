#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// filter Mutect2 calls with GATK
process FilterMutectCalls {
    publishDir "${params.outdir}/filtered", mode: 'copy'

    input:
    path unfiltered_vcf
    path mutect_idx
    path mutect_idx_fai
    path mutect_dict
    path vcf_stats

    output:
    file("${unfiltered_vcf.baseName}_filtered.vcf")

    script:
    """
    gatk FilterMutectCalls \
    -R ${mutect_idx} \
    -V ${unfiltered_vcf} \
    -O ${unfiltered_vcf.baseName}_filtered.vcf \
    --read-index ${mutect_idx_fai} \
    --sequence-dictionary ${mutect_dict} \
    --stats ${vcf_stats}
    """
}