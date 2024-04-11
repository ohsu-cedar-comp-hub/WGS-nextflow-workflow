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
    path read_orientation_model
    path segmentation_tabel
    path contamination_table

    output:
    file("${unfiltered_vcf.baseName}_filtered.vcf")

    script:
    """
    gatk FilterMutectCalls \
    -R ${mutect_idx} \
    -V ${unfiltered_vcf} \
    --tumor-segmentation ${params.segmentation_table} \
    --contamination-table ${params.contamination_table} \
    -O ${unfiltered_vcf.baseName}_filtered.vcf \
    --read-index ${mutect_idx_fai} \
    --sequence-dictionary ${mutect_dict} \
    --ob-priors ${params.read_orientation_model} \
    --stats ${vcf_stats}
    """

}

workflow {
    // define input parameters
    unfiltered_vcf = file(params.unfiltered_vcf)
    mutect_idx = file(params.mutect_idx)
    mutect_idx_fai = file(params.mutect_idx_fai)
    mutect_dict = file(params.mutect_dict)
    ob_priors = file(params.read_orientation_model)
    vcf_stats = file(params.vcf_stats)

    // run the FilterMutectCalls process
    FilterMutectCalls(unfiltered_vcf, mutect_idx, mutect_idx_fai, mutect_dict, vcf_stats)
}

