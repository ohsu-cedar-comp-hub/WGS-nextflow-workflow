
process CGPWGS {
    // execute the container exactly how it is set up in the sanger evotypes pipeline : ascat, caveman, brass included

    input:

    output:

    script:
    """
    ds-cgpwgs.pl \
    -c $NSLOTS \
    -r /var/spool/ref/core_ref_GRCh38_hla_decoy_ebv.tar.gz \
    -a /var/spool/ref/VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91.tar.gz \
    -si /var/spool/ref/SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment.tar.gz \
    -cs /var/spool/ref/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz \
    -qc /var/spool/ref/qcGenotype_GRCh38_hla_decoy_ebv.tar.gz \
    -e 'chrUn%,HLA%,%_alt,%_random,chrM,chrEBV' \
    -t /var/spool/tmdata/${tumorx}.bam \
    -tidx /var/spool/tmdata/${tumorx}.bam.bai \
    -n /var/spool/nmdata/${normalx}.bam \
    -nidx /var/spool/nmdata/${normalx}.bam.bai \
    -o /var/spool/results
    """
}



