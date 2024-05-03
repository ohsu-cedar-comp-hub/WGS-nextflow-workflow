#!/usr/bin/env nextflow

include {MutectSplitByChromosome} from '../tools/mutect/mutect_split_by_chromosome.nf'

workflow {
    // Run Mutect2
    MutectSplitByChromosome(file(params.tumor_bam) ,file(params.normal_bam), file(params.bed_files).collect(), val(params.id), file(params.mutect_idx), file(params.pon), "test")
}

