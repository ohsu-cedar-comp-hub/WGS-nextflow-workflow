#!/usr/bin/env nextflow

// filter Mutect2 calls with GATK
process FILTERMUTECT {
    // Set maximum memory
    // memory '40 GB'

    container "${params.container_gatk}"

    publishDir "${params.outdir}/filtered", mode: 'copy'

    input:
    path unfiltered_vcf
    path unfiltered_vcf_index
    path mutect_idx
    // path mutect_idx_fai
    // path mutect_dict
    path vcf_stats
    // path read_orientation_model
    path segmentation_table
    path contamination_table
    val sample_id

    output:
    path("${sample_id}_filtered.vcf")

    script:
    """
    gatk FilterMutectCalls \
    -R ${params.mutect_idx} \
    -V ${unfiltered_vcf} \
    --tumor-segmentation ${segmentation_table} \
    --contamination-table ${contamination_table} \
    -O ${sample_id}_filtered.vcf \
    --stats ${vcf_stats}
    """
// --ob-priors ${read_orientation_model} \    
// --read-index ${mutect_idx_fai} \
// --sequence-dictionary ${mutect_dict} \
}
