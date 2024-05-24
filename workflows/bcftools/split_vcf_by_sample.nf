#!/usr/bin/env nextflow

include {SplitVcfBySample} from '../../tools/bcftools/split_vcf_by_sample.nf'

// Define the workflow
workflow {
    // Run SplitVcfBySample
       SplitVcfBySample(file(params.annotated_vcf), "test")
}
