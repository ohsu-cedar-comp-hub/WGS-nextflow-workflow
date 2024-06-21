#!/usr/bin/env nextflow

files_ch = Channel.fromPath("${params.bucket_dir}*.${params.filesuffix}")

include { TOFASTQ } from '../../tools/samtools/to_fastq.nf'

workflow {
    TOFASTQ(files_ch)
    files_ch.view()
}