#!/usr/bin/env nextflow

// Create Channel 
// all_pairs_ch is a read pairs channel structured like [id, [r1.fq, r2.fq]]
all_pairs_ch = Channel.fromFilePairs(params.all_read_pairs)

// import modules 
include { FASTQC as FASTQCRAW} from '../../tools/qc/fastqc/fastqc.nf'
include { FASTQC as FASTQCTRIM } from '../../tools/qc/fastqc/fastqc.nf'
include { TRIMMOMATICPE } from '../../tools/trimmomatic/trimmomatic.nf'
include { MULTIQC } from '../../tools/qc/multiqc/multiqc.nf'
include { BWAMEM2 } from '../../tools/bwa/bwamem2.nf'
include { SORT; SORTANDINDEX } from '../../tools/samtools/sort_and_index.nf'
include { MARKDUPLICATES } from '../../tools/gatk/mark_duplicates.nf'

workflow 
workflow {
    // trimmomatic
    TRIMMOMATICPE(all_pairs_ch, params.truseq3pefile, params.outdir)

    // fastqc on raw and trimmed reads 
    FASTQCRAW(all_pairs_ch, params.outdir)
    FASTQCTRIM(TRIMMOMATICPE.out.trim_reads, params.outdir)

    // gather all files output from fastqc processes
    multi_ch = FASTQCRAW.out.zip.mix(FASTQCTRIM.out.zip).collect()
    
    // pass to multiqc
    MULTIQC(multi_ch, all_pairs_ch)

    // align with bwa-mem2
    BWAMEM2(TRIMMOMATICPE.out.trim_reads, params.idx)

    // sort with samtools 
    SORT(BWAMEM2.out)

    // mark duplicates
    MARKDUPLICATES(SORT.out)

    // sort and index with samtools to prep for gatk somatic variant calling
    SORTANDINDEX(MARKDUPLICATES.out.bam)
}
