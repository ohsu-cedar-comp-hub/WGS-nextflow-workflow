process LEARNORIENTATION {
    publishDir "${params.outdir}/svc/f1r2files/", mode: 'copy'
    
    container "${params.container_gatk}"

    input:
    path f1r2file
    val sample_id

    output:
    path("${sample_id}_read-orientation-model.tar.gz")

    script:
    """
    gatk LearnReadOrientationModel -I ${f1r2file.join(' -I ')} -O ${sample_id}_read-orientation-model.tar.gz
    """
}
