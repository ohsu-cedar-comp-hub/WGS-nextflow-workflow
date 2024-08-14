#!/usr/bin/env nextflow

/*
-o Folder to output result to.
-t Tumour BAM/CRAM file
-n Normal BAM/CRAM file
-sg Snp GC correction file
-r Reference fasta
-q Minimum base quality required before allele is used. Default 20
-g Sample gender (XX, XY, L). When 'L' see '-l'
-l Using a list of loci, default when '-L' [share/gender/GRCh37d5_Y.loci]. These are loci that will not present at all in a female sample
-rs Reference species [BAM HEADER]
-ra Reference assembly [BAM HEADER]
-pr Sequencing protocol (e.g. WGS, WXS)
-pl Seqeuncing platform [BAM HEADER]
-c Number of cores to use. Default 1. Recommend max 2 during 'input' process.
-force Force completion - solution not possible. 
    Adding this will result in successful completion of analysis even when ASCAT can't generate a solution.  
    A default copynumber of 5/2 (tumour/normal) and contamination of 30% will be set along with a comment 
    in '*.samplestatistics.csv' to indicate this has occurred.

OPTIONAL:
-pu tumor purity and -pi tumor ploidy, not included in cgpwgs_NOORANI_corrected_111-116.sh https://github.com/elisabethgoldman/sanger-evotypes/blob/main/cgpwgs_NOORANI_corrected_111-116.sh
*/

process ASCAT {
    container "${params.container_cgpwgs}"
    publishDir "${params.outdir}/cgp-wgs/ascat", mode: 'copy'

    cpus = 1

    input:
    path tumor_bam
    path normal_bam
    path snpgccorrections_tsv
    path reference_fa
    path gender_tsv

    output:
    path ("something tumor cn"), emit: tumor_cn
    path ("something normal cn"), emit: normal_cn

    script:
    """
    ascat.pl \
    -o . \
    -t ${tumor_bam} \
    -n ${normal_bam} \
    -sg ${snpgccorrections_tsv} \
    -r ${reference_fa} \
    -q 20 \
    -g L \
    -l ${gender_tsv} \
    -rs HUMAN \
    -ra GRCh38.d1.vd1 \
    -pr WGS \
    -pl illumina \
    -c 1 \
    -force
    """
}

workflow {
    ASCAT(tumor_bam, normal_bam, params.snpgccorrections_tsv, params.reference_fa, params.gender_tsv)
}