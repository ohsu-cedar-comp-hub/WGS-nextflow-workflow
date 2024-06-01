process LEARNORIENTATION {
    publishDir "${params.outdir}/svc/f1r2files/", mode: 'copy'
    
    input:
    val f1r2file
    val sample_id

    output:
    path 'read-orientation-model.tar.gz'

    script:
    """
    gatk LearnReadOrientationModel ${f1r2file.join('-I ')} -O read-orientation-model.tar.gz
    """
}
