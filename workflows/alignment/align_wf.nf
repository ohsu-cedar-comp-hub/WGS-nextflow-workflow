#!/usr/bin/env nextflow

// Create Channel 
all_pairs_ch = Channel.fromFilePairs(params.all_read_pairs)

// import modules 
include { FASTQC as FASTQCRAW} from '../../tools/qc/fastqc/fastqc.nf'
include { FASTQC as FASTQCTRIM } from '../../tools/qc/fastqc/fastqc.nf'
include { TRIMMOMATICPE } from '../../tools/trimmomatic/trimmomatic.nf'
include { MULTIQC } from '../../tools/qc/multiqc/multiqc.nf'
include { BWAMEM2 } from '../../tools/bwa/bwamem2.nf'

// workflow 
workflow {
    // all_pairs_ch is a read pairs channel structured like [id, [r1.fq, r2.fq]]

    // trimmomatic
    TRIMMOMATICPE(all_pairs_ch, params.truseq3pefile, params.outdir)

    // fastqc on raw and trimmed reads 
    FASTQCRAW(all_pairs_ch, params.outdir)
    FASTQCTRIM(TRIMMOMATICPE.out.trim_reads, params.outdir)

    // create a channel that is a flat list of all the fastqc zip files
    multi_ch = FASTQCRAW.out.zip.mix(FASTQCTRIM.out.zip).flatten()
    // pass to multiqc
    MULTIQC(multi_ch)

    // align with bwa-mem2
    BWAMEM2(TRIMMOMATICPE.out.trim_reads, params.idx, params.id)
}
