#!/usr/bin/env nextflow

// filter Mutect2 calls with GATK
process FilterMutectCalls {
    // Set maximum memory
    memory '40 GB'

    publishDir "${params.outdir}/vcf_filtered", mode: 'copy'

    input:
    path mutect_idx
    path unfiltered_vcf
    path read_orientation_model
    path contamination_table
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
    --ob-priors ${read_orientation_model} \
    --contamination-table ${contamination_table} \
    --read-index ${mutect_idx_fai} \
    --sequence-dictionary ${mutect_dict} \
    --stats ${vcf_stats}
    """

}
