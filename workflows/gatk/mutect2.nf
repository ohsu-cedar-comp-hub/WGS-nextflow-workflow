#!/usr/bin/env nextflow

include {Mutect2} from '../../tools/gatk/mutect2.nf'

workflow {
    // Run Mutect2
    Mutect2(file(params.tumor_bam_sorted) ,file(params.normal_bam_sorted), file(params.mutect_idx), file(params.pon), file(params.gnomad), params.ID)
}
