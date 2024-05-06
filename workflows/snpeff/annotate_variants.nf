#!/usr/bin/env nextflow

include {Annotate_Variants} from '../../tools/snpeff/annotate_variants.nf'

// Define the workflow for annotating variants
workflow {

    Annotate_Variants (file(params.filtered_vcf), "test")

}
