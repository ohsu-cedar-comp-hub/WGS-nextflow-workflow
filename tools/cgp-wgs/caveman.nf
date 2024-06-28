#!/usr/bin/env nextflow

process CAVEMAN {
    container "${params.container_cgpwgs}"
    publishDir "${params.outdir}/cgp-wgs/caveman", mode: 'copy'

    input:
    ascat bed copy number output tumor
    ascat bed copy number output normal
    tumor bam
    normal bam
    all reference files

    output:
    vcf 

    script:
    """
    export REF=/referenceArea
    export POUT=/output
    mkdir $POUT/result
    export CAVE=/exampleData

    caveman.pl \
    -reference $REF/genome.fa.fai \
    -outdir $POUT/result \
    -tumour-bam $CAVE/tumour/HCC1143.bam \
    -normal -bam $CAVE/normal/HCC1143_BL.bam \
    -ignore -file $REF/ignore_regions.bed \
    -tumour-cn $CAVE/tumour/HCC1143.cn.bed \
    -normal-cn $CAVE/normal/HCC1143_BL.cn.bed \
    -species Human \
    -species-assembly GRCh37d5 \
    -flag-bed-files $REF/flagging \
    -germline-indel $CAVE/tumour/HCC1143_vs_HCC1143_BL.germline.bed \
    -unmatched-vcf $REF/flagging/unmatched_vcf \
    -seqType genomic \
    -normal-contamination \
    -threads 32 \
    -limit 32 \
    -flagConfig $REF/flagging/ \
    -flagToVcfConfig $REF/flagging/ \
    -annot-bed-files $REF/vagrent/ \
  >& $POUT/caveman_run.log
    """
}