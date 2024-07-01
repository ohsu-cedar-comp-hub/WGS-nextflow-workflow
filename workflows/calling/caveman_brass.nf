#!/usr/bin/env nextflow

// create channels


workflow {
    ascat
    caveman(ascat.out, reference files)
    pindel()
    brass()
}