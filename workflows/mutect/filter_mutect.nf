#!/usr/bin/env nextflow

include {FilterMutectCalls} from '../tools/mutect/filter_mutect.nf

// define workflow
workflow {
    // run the FilterMutectCalls process
    FilterMutectCalls(file(params.mutect_idx),file(params.unfiltered_vcf), file(params.read_orientation_model), file(params.contamination_table), file(params.mutect_idx_fai), file(params.mutect_dict), file(params.vcf_stats), "test")
}


