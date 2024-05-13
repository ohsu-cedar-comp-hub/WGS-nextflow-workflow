#!/usr/bin/env nextflow

include {SiftVariants} from '../../tools/snpeff/sift_variants.nf'

// Define the workflow for annotating variants
workflow {

    SiftVariants (file(params.filtered_vcf), params.ID)

}