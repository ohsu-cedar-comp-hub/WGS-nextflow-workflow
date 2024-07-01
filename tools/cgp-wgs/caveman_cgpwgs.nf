#!/usr/bin/env nextflow

/* caveman.pl arguments:

-r Path to genome.fa.fai file and associated .fa file.
-ig Path to ignore file. 1-based first coordinate bed style format of regions for caveman not to analyse.
-b Bed file location for flagging (eg dbSNP.bed NB must be sorted.)
-ab Annotation BED files - required for pulldown/WXS
-u Directory containing unmatched normal VCF files
-s Species name for (output in VCF) e.g HUMAN
-sa Species assembly for (output in VCF) e.g. 37
-t Number of threads allowed on this machine
-st Sequencing type (pulldown|exome|genome|genomic|followup|targeted|rna_seq) - Passed to flagging
-tc Path to tumour copy number file (1-based first coordinate bed style format). All analysed bases must have a CN assigned.
-nc Path to normal copy number file (1-based first coordinate bed style format). All analysed bases must have a CN assigned.
-td 5 Use empty file with -td for blanket copynumber
-nd 2 Use empty file with -nd for blanket copynumber
-tb Path to mapped, indexed, duplicate marked/removed tumour bam file.
-nb Path to mapped, indexed, duplicate marked/removed normal bam file.
-c Config ini file to use for flag list and settings
-f Config::Inifiles style config file containing VCF flag code to flag name conversions
-e Modify the read count threshold used when splitting into mstep/estep jobs (from ds-cgpwgs.pl: Reads target per caveman split section. Default 350000.)
-o Directory to write output to. During processing a temp folder will be generated in this area, should the process fail only delete this if you are unable to resume the process.
Final output files are: muts.vcf.gz, snps.vcf.gz, no_analysis.bed.gz no_analysis.bed.gz.tbi
-x contig exclude (lacking documentation)
-p Only process this step then exit
*/

process CAVEMAN_SETUP {
    container "${params.container_cgpwgs}"
    publishDir "${params.outdir}/cgp-wgs/caveman", mode: 'copy'

    cpu = 4

    input:
    path caveman_ref_fai
    path ignore_hi_depth_tsv
    path tumor_cn
    path normal_cn
    path tumor_bam
    path normal_bam

    output:
    path ("caveman config caveman.cfg.ini"), emit: config
    path ("alg_bean something"), emit: alg_bean

    // script from analysisWGS.sh at https://github.com/cancerit/dockstore-cgpwgs/blob/develop/scripts/analysisWGS.sh
    script:
    """
    caveman.pl \
    -r ${caveman_ref_fai} \
    -ig ${ignore_hi_depth_tsv} \
    -b $REF_BASE/caveman/flagging \ ***?????????***
    -ab $REF_BASE/vagrent \ ***?????????***
    -u $REF_BASE/caveman \ ***?????????***
    -s HUMAN \
    -sa GRCh38.d1.vd1 \
    -t 4 \
    -st WGS \
    -tc ${tumor_cn} \ 
    -nc ${normal_cn} \
    -td 5 -nd 2 \ ***?????????***
    -tb ${tumor_bam} \
    -nb ${normal_bam} \
    -c flag.vcf.config.WGS.ini \ ***?????????***
    -f flag.to.vcf.convert.ini \ ***?????????***
    -e 350000 \
    -o . \
    -x $CONTIG_EXCLUDE \ ***?????????***
    -p setup
    """
}

process CAVEMAN_SPLIT {
    container "${params.container_cgpwgs}"
    publishDir "${params.outdir}/cgp-wgs/caveman", mode: 'copy'

    cpu = 4 ??

    input:
    path caveman_ref_fai
    path ignore_hi_depth_tsv
    path tumor_cn
    path normal_cn
    path tumor_bam
    path normal_bam

    output:


    script:
    """
    caveman.pl \
    -r $REF_BASE/genome.fa.fai \
    -ig $REF_BASE/caveman/HiDepth.tsv \
    -b $REF_BASE/caveman/flagging \
    -ab $REF_BASE/vagrent \
    -u $REF_BASE/caveman \
    -s '$SPECIES' \
    -sa $ASSEMBLY \
    -t $CPU \
    -st $PROTOCOL \
    -tc $TMP/tum.cn.bed \
    -nc $TMP/norm.cn.bed \
    -td 5 -nd 2 \
    -tb $BAM_MT_TMP \
    -nb $BAM_WT_TMP \
    -c $SNVFLAG \
    -f $REF_BASE/caveman/flagging/flag.to.vcf.convert.ini \
    -e $CAVESPLIT \
    -o $OUTPUT_DIR/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman \
    -x $CONTIG_EXCLUDE \
    -p split
    """
}

process CAVEMAN_MAIN {

    cpus: Recommended 32 cores (minimum)
    memory: 4GB per core

    script:
    """
    caveman.pl \
 -r $REF_BASE/genome.fa.fai \
 -ig $REF_BASE/caveman/HiDepth.tsv \
 -b $REF_BASE/caveman/flagging \
 -ab $REF_BASE/vagrent \
 -u $REF_BASE/caveman \
 -s '$SPECIES' \
 -sa $ASSEMBLY \
 -t $CPU \
 -st $PROTOCOL \
 -tc $TMP/tum.cn.bed \
 -nc $TMP/norm.cn.bed \
 -td 5 -nd 2 \
 -tb $BAM_MT_TMP \
 -nb $BAM_WT_TMP \
 -c $SNVFLAG \
 -f $REF_BASE/caveman/flagging/flag.to.vcf.convert.ini \
 -e $CAVESPLIT \
 -o $OUTPUT_DIR/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman \
 -x $CONTIG_EXCLUDE \
 -k $NORM_CONTAM \
 -no-flagging -noclean"
    """
}

process CAVEMAN_FLAG {

    script:

    """
    caveman.pl \
 -r $REF_BASE/genome.fa.fai \
 -ig $REF_BASE/caveman/HiDepth.tsv \
 -b $REF_BASE/caveman/flagging \
 -ab $REF_BASE/vagrent \
 -u $REF_BASE/caveman \
 -s '$SPECIES' \
 -sa $ASSEMBLY \
 -t $CPU \
 -st $PROTOCOL \
 -tc $TMP/tum.cn.bed \
 -nc $TMP/norm.cn.bed \
 -td 5 -nd 2 \
 -tb $BAM_MT_TMP \
 -nb $BAM_WT_TMP \
 -c $SNVFLAG \
 -f $REF_BASE/caveman/flagging/flag.to.vcf.convert.ini \
 -e $CAVESPLIT \
 -o $OUTPUT_DIR/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman \
 -x $CONTIG_EXCLUDE \
 -k $NORM_CONTAM \
 -in $GERMLINE_BED.gz \
 -p flag"
    """
}

workflow {
    CAVEMAN_SETUP(caveman_ref_fai, ignore_hi_depth_tsv, ASCAT.out.tumor_cn, ASCAT.out.normal_cn, tumor_bam, normal_bam, )
    CAVEMAN_SPLIT()
    CAVEMAN_MAIN()
    CAVEMAN_FLAG()
}