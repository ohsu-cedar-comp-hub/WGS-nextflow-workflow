#!/usr/bin/env nextflow

include {mutect2} from '../tools/mutect/mutect2.nf'

workflow {
    // Run Mutect2
    mutect2(file(params.tumor_bam) ,file(params.normal_bam), file(params.idx), file(params.pon), file(params.gnomad), params.ID)
}
