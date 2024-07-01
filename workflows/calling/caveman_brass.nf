#!/usr/bin/env nextflow

// create channels


workflow {
    ASCAT(tumor_bam, normal_bam, params.snpgccorrections_tsv, params.reference_fa, params.gender_tsv)
    caveman(ascat.out, reference files)
    pindel()
    brass()
}