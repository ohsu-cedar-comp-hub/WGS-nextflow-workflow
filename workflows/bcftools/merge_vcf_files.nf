#!/usr/bin/env nextflow

include {MergeFiles} from '../../tools/bcftools/merge_vcf_files.nf


workflow {
   MergeFiles(file(params.outdir/svc).collect(), val(params.id))


}
