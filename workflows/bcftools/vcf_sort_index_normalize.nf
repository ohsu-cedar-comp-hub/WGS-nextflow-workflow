#!/usr/bin/env nextflow

include {VcfSortIndexNormalize} from '../../tools/bcftools/vcf_sort_index_normalize.nf'

// Define the workflow
workflow {
    // Run vcfSortIndexNormalize 
    VcfSortIndexNormalize(file(params.unfiltered_vcf), params.ID)
}
