#!/usr/bin/env nextflow

// define workflow
workflow {
    // define input parameters
    unfiltered_vcf = file(params.unfiltered_vcf)
    mutect_idx = file(params.mutect_idx)
    mutect_idx_fai = file(params.mutect_idx_fai)
    mutect_dict = file(params.mutect_dict)
    vcf_stats = file(params.vcf_stats)

    // run the FilterMutectCalls process
    FilterMutectCalls(unfiltered_vcf, mutect_idx, mutect_idx_fai, mutect_dict, vcf_stats)
}
