#!/usr/bin/env nextflow

// Process for learning gatk read orientation model
process LearnReadOrientationModel {
    // Set maximum memory
    memory '40 GB'

    publishDir "${params.outdir}/tables", mode: 'copy'

    input: // if multiple tumor samples are being processed, only a single f1r2_tar_gz output is necessary, which will contain info for all included samples
    path f1r2_tar_gz
    val ID

    output:
    path "read_orientation_model.tar.gz"

    script:
    """
    gatk LearnReadOrientationModel -I ${f1r2_tar_gz} -O read_orientation_model.tar.gz
    """
}
