#!/usr/bin/env nextflow

// Parameters 
params.normal_reads = "/home/groups/CEDAR/lancasru/WGS_COH_NF/nextflow_test/references/sliced_fastqs/*_G_R{1,2}.fastq.gz"
params.tumor_reads = "/home/groups/CEDAR/lancasru/WGS_COH_NF/nextflow_test/references/sliced_fastqs/*_T_R{1,2}.fastq.gz"
params.all_reads = "/home/groups/CEDAR/lancasru/WGS_COH_NF/nextflow_test/references/sliced_fastqs/*.fastq.gz"
params.all_read_pairs = "/home/groups/CEDAR/lancasru/WGS_COH_NF/nextflow_test/references/sliced_fastqs/*_R{1,2}.fastq.gz"
params.outdir = "/home/groups/CEDAR/lancasru/WGS_COH_NF/nextflow_test/channel_snippet"
params.truseq3pefile = "/home/groups/CEDAR/lancasru/WGS_COH_NF/nextflow_test/references/TruSeq3-PE.fa"

// Create Channels 

// Used for FASTQC
all_ch = Channel.fromPath(params.all_reads)

// Used for TRIMMOMATICPE
all_pairs_ch = Channel.fromFilePairs(params.all_read_pairs)

normal_ch = Channel.fromFilePairs(params.normal_reads)
tumor_ch = Channel.fromFilePairs(params.tumor_reads)

// import modules 

include { FASTQC } from '../../tools/qc/fastqc/fastqc.nf'
include { TRIMMOMATICPE } from '../../tools/trimmomatic/trimmomatic.nf'
// workflow 


workflow {
    FASTQC(all_ch, params.outdir)
    TRIMMOMATICPE(all_pairs_ch, params.truseq3pefile, params.outdir)
}
