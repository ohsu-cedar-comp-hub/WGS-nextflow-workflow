bam_file = Channel.fromPath("/home/exacloud/gscratch/CEDAR/grieco/nextflow/H_51720/aligned/*.bam")
flat_ch= bam_file.flatten()

// import modules 
include { SORT; SORTANDINDEX } from '../../tools/samtools/sort_and_index.nf'
include { MARKDUPLICATES } from '../../tools/gatk/mark_duplicates.nf'

workflow {

    // sort with samtools 
    SORT(flat_ch)

    // mark duplicates
    MARKDUPLICATES(SORT.out)
    
    // sort and index with samtools to prep for gatk somatic variant calling
    SORTANDINDEX(MARKDUPLICATES.out.bam)
}