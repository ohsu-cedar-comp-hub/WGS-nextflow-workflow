#!/usr/bin/env nextflow

workflow {
    // define input paramaters for Mutect2
    tumor_bam=file(params.tumor_bam)
    normal_bam=file(params.normal_bam)
    idx=file(params.mutect_idx)
    pon=file(params.pon)
    // Run Mutect2
    mutect2(tumor_bam, normal_bam, idx, pon)
}
