#!/usr/bin/env nextflow

// QC, trimmomatic (optional), alignment with bwa-mem2, and mark duplicates

// import tool modules

include {FastQC} from '../../tools/qc/fastqc/fastqc.nf'
include {TrimmomaticPE} from '../../tools/trimmomatic/trimmomatic.nf'
include {MultiQC} from '../../tools/qc/multiqc/multiqc.nf'
include {bwaMem2Alignment} from '../../tools/bwa/bwamem2.nf'
include {SortAndIndex} from '../../tools/samtools/sort_and_index.nf'
include {MarkDuplicates} from '../../tools/gatk/mark_duplicates.nf'
include {SortMarkedDuplicates} from '../../tools/samtools/sort_marked_duplicates.nf'

// define the workflow

workflow {
    FastQC(params.read1, params.read2, params.ID)
    TrimmomaticPE(params.read1, params.read2, params.truseq3pefile, params.ID)
    FastQC(file(params.trim_read1), file(params.trim_read2), params.ID)
    MultiQC(file(params.fastqc_read1), file(params.fastqc_read2), file(params.trim_fastqc_read1), file(params.trim_fastqc_read2), params.ID)
    bwaMem2Alignment(params.trim_read1, params.trim_read2, params.idx, params.ID)
    SortAndIndex(params.bam_unsorted, params.ID)
    MarkDuplicates(params.tumor_bam_sorted, params.normal_bam_sorted, params.ID)
    SortMarkedDuplicates (params.bam_duplicates_unsorted, params.ID)
}

