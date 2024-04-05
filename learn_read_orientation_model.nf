// Process for learning gatk read orientation model
process LearnReadOrientationModel {
    publishDir "${params.outdir}/models", mode: 'copy'

    input:
    path f1r2_tar_gz

    output:
    path "read-orientation-model.tar.gz"

    script:
    """
    gatk LearnReadOrientationModel -I ${f1r2_tar_gz} -O read_orientation_model.tar.gz
    """
}
// define workflow
workflow {
    // Define input parameters
    f1r2_tar_gz = file(params.f1r2_tar_gz)

    // run process
    LearnReadOrientationModel(f1r2_tar_gz)
}
