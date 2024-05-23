#!/usr/bin/env nextflow

// Parameters 
params.normal_reads = "/home/groups/CEDAR/lancasru/WGS_COH_NF/nextflow_test/references/sliced_fastqs/*_G_R{1,2}.fastq.gz"
params.tumor_reads = "/home/groups/CEDAR/lancasru/WGS_COH_NF/nextflow_test/references/sliced_fastqs/*_T_R{1,2}.fastq.gz"
params.all_reads = "/home/groups/CEDAR/lancasru/WGS_COH_NF/nextflow_test/references/sliced_fastqs/*.fastq.gz"
params.all_read_pairs = "/home/groups/CEDAR/lancasru/WGS_COH_NF/nextflow_test/references/sliced_fastqs/*_R{1,2}.fastq.gz"
params.outdir = "/home/groups/CEDAR/lancasru/WGS_COH_NF/nextflow_test/channel_snippet"

// Channel creation 

// Used for FASTQC
all_ch = Channel.fromPath(params.all_reads)

normal_ch = Channel.fromFilePairs(params.normal_reads)
tumor_ch = Channel.fromFilePairs(params.tumor_reads)
all_pairs_ch = Channel.fromFilePairs(params.all_read_pairs)


process FASTQC {
    // debug true
    
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    path reads
    path outdir

    output:
    path("*.zip"), emit: zip
    path("*.html"), emit: html
    
    script:
    """
    if [ -d $outdir/fastqc ]; then
        :
    else
        mkdir $outdir/fastqc
    fi 
    
    /usr/local/FastQC/fastqc $reads
    """
}

workflow {
    FASTQC(all_ch, params.outdir)
}

/*
process TRIMMOMATIC {
    debug true 

    input:
    path read 1
    path read 2
    path truseq3pefile
    path outdir

    output:

    script: 
    """

    """


}
*/

// trim_ch = Channel.fromFilePairs(params.all_read_pairs)
// use mapping for this to create two channels; one for read 1 and one for read 2 and ensure they are fed in as queue channel with tag so they dont get mixed up 

/*
workflow {
    FASTQC(all_ch, params.outdir)
    TRIMMOMATIC(trim_ch, params.outdir, params.truseq3pefile)
}
*/
