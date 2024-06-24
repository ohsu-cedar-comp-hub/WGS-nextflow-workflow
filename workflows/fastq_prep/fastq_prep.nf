#!/usr/bin/env nextflow

// before running, make sure the temp directory is set in the same dir you launch the sbatch from, otherwise samtools collate cannot locate the correct temp dir
// > cd <work directory>
// > export TMPDIR=`pwd`/tmp

files_ch = Channel.fromPath("${params.bucket_dir}/*.${params.filesuffix}")

include { TOFASTQSORT; TOFASTQ } from '../../tools/samtools/to_fastq.nf'

workflow {
    TOFASTQSORT(files_ch)
    sorted = TOFASTQSORT.out.collect()
    TOFASTQ(sorted)
}