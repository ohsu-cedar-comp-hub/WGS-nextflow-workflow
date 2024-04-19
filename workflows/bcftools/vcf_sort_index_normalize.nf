#!/usr/bin/env nextflow

include {vcfSortIndexNormalize} from '../../tools/bcftools/vcf_sort_index_normalize.nf'

// Define the workflow
workflow {
    // Run vcfSortIndexNormalize 
    bwaMem2Alignment(file(params.unfiltered_vcf), "test")
}
