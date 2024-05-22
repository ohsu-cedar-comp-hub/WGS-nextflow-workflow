#!/usr/bin/env nextflow

// QC, trimmomatic (optional), alignment with bwa-mem2, and mark duplicates

// import tool modules

include {FastQC as FastQCRaw} from '../../tools/qc/fastqc/fastqc.nf'
include {FastQC as FastQCTrim} from '../../tools/qc/fastqc/trim_fastqc.nf'
include {TrimmomaticPE} from '../../tools/trimmomatic/trimmomatic.nf'
include {MultiQC} from '../../tools/qc/multiqc/multiqc.nf'
include {BwaMem2Alignment} from '../../tools/bwa/bwamem2.nf'
include {SortAndIndex} from '../../tools/samtools/sort_and_index.nf'
include {MarkDuplicates} from '../../tools/gatk/mark_duplicates.nf'
include {SortMarkedDuplicates} from '../../tools/samtools/sort_marked_duplicates.nf'

// define the workflow

workflow {
    FastQCRaw(params.read1, params.read2, params.ID)
    TrimmomaticPE(params.read1, params.read2, params.truseq3pefile, params.ID)
    FastQCTrim(params.trim_read1, params.trim_read2, params.ID)
    MultiQC(params.fastqc_read1, params.fastqc_read2, params.trim_fastqc_read1, params.trim_fastqc_read2, params.ID)
    BwaMem2Alignment(params.trim_read1, params.trim_read2, params.idx, params.ID)
    SortAndIndex(params.bam_unsorted, params.ID)
    MarkDuplicates(params.bam_sorted, params.ID)
    SortMarkedDuplicates(params.bam_duplicates_unsorted, params.ID)
}

