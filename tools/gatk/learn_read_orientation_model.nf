#!/usr/bin/env nextflow

// Process for learning gatk read orientation model
process LearnReadOrientationModel {
    publishDir "${params.outdir}/tables", mode: 'copy'

    input: // if multiple tumor samples are being processed, only a single f1r2_tar_gz output is necessary, which will contain info for all included samples
    path f1r2_tar_gz

    output:
    path "read_orientation_model.tar.gz"

    script:
    """
    gatk LearnReadOrientationModel -I ${params.f1r2_tar_gz} -O read_orientation_model.tar.gz
    """
}
// define workflow
workflow {
    // Define input parameters
    f1r2_tar_gz = file(params.f1r2_tar_gz)

    // run process
    LearnReadOrientationModel(f1r2_tar_gz)
}
