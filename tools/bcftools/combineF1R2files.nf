process ReadOrientationModel {
    publishDir "${params.outdir}/svc/f1r2files/", mode: 'copy'
    
    input:
    val chromosomes
    path f1r2file

    output:
    path 'read-orientation-model.tar.gz'

    script:
    """

    all_f1r2_input=""
    for chromosome in ${chromosomes}; do
        all_f1r2_input+=" -I \${chromosome}.f1r2.tar.gz"
    done

    LearnReadOrientationModel \$all_f1r2_input -O read-orientation-model.tar.gz
    """
}
