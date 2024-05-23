#!/usr/bin/env nextflow

// Parameters for params file 
params.all_read_pairs = "/home/groups/CEDAR/lancasru/WGS_COH_NF/nextflow_test/references/sliced_fastqs/*_R{1,2}.fastq.gz"
params.outdir = "/home/groups/CEDAR/lancasru/WGS_COH_NF/nextflow_test/channel_snippet"
params.truseq3pefile = "/home/groups/CEDAR/lancasru/WGS_COH_NF/nextflow_test/references/TruSeq3-PE.fa"
params.idx = "/home/groups/CEDAR/goldmael/projects/align/resources/bwamem2/GRCh38.bm2.d1"
params.id = "UUID-UUID-UUID"

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
