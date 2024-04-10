#!/usr/bin/env nextflow

workflow {
    // define input paramaters for Mutect2
    input_normal_file=file(params.normal_bam)
    input_tumor_file=file(params.tumor_bam)
    idx=file(params.mutect_idx)
    pon=file(params.pon)
    // Run Mutect2
    mutect2(input_normal_file, input_tumor_file, idx, pon)
}
