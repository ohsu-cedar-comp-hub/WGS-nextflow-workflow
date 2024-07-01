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
tumor purity and tumor ploidy, not included in cgpwgs_NOORANI_corrected_111-116.sh https://github.com/elisabethgoldman/sanger-evotypes/blob/main/cgpwgs_NOORANI_corrected_111-116.sh
*/

process ASCAT {
    container "${params.container_cgpwgs}"
    publishDir "${params.outdir}/cgp-wgs/ascat", mode: 'copy'

    input:
    tumor bam
    normal bam

    output:
    tumor copy number bed
    normal copy number bed

    script:
    """
    ascat.pl \
    -o $OUTPUT_DIR/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/ascat \
    -t $BAM_MT_TMP \
    -n $BAM_WT_TMP \
    -sg $REF_BASE/ascat/SnpGcCorrections.tsv \
    -r $REF_BASE/genome.fa \
    -q 20 \
    -g L \
    -l $REF_BASE/gender.tsv \
    -rs '$SPECIES' \
    -ra $ASSEMBLY \
    -pr $PROTOCOL \
    -pl ILLUMINA \
    -c $CPU \
    -force
    """
}