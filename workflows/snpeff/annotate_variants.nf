#!/usr/bin/env nextflow

include {AnnotateVariants} from '../../tools/snpeff/annotate_variants.nf'

// Define the workflow for annotating variants
workflow {

    AnnotateVariants (file(params.filtered_vcf), params.ID)

}
