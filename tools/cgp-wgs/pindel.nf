#!/usr/bin/env nextflow


process PINDEL {
    script:
    """
    do_parallel[cgpPindel]="pindel.pl \
 -o $OUTPUT_DIR/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/pindel \
 -r $REF_BASE/genome.fa \
 -t $BAM_MT_TMP \
 -n $BAM_WT_TMP \
 -s $REF_BASE/pindel/simpleRepeats.bed.gz \
 -u $REF_BASE/pindel/pindel_np.gff3.gz \
 -f $REF_BASE/pindel/${PROTOCOL}_Rules.lst \
 -g $REF_BASE/vagrent/codingexon_regions.indel.bed.gz \
 -st $PROTOCOL \
 -as $ASSEMBLY \
 -sp '$SPECIES' \
 -e $CONTIG_EXCLUDE \
 -b $REF_BASE/pindel/HiDepth.bed.gz \
 -c $PINDEL_CPU \
 -sf $REF_BASE/pindel/softRules.lst"
    """
}