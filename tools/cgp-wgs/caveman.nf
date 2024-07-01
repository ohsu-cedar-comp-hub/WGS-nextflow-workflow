#!/usr/bin/env nextflow

// utilize scatter-gather for optimized parallel processing 
// details here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6097605/#S12title


// Run the single job setup step
process CAVEMAN_SETUP {

}

// Run the split step once per contig in the *.fa.fai file:
process CAVEMAN_SPLIT {


    script:
    """
    ...

    -process split \
    –index ${contig_channel; 1,2,3...n. run once per contig in fa.fai file (eg, if 50 lines long, 50 jobs spawn)}
    """
}

// Run a single merge split sections job to create a single splitList file:
process CAVEMAN_CONCAT {
    script:
    """
    ...
    -process split_concat \
    –index 1
    """
}

// run a mstep job per section in the splitList file
process CAVEMAN_MSTEP {
    script:
    """
    ...
    -process mstep \
    –index ${splitlist_channel; 1,2,3...n. Calculate the number of sections in the split file and run a mstep job per section}
    """
}

// Merge the profiles generated in the mstep jobs into a single file
process CAVEMAN_MERGE {
    script:
    """
    ...
    -process merge \
    –index 1
    """
}

// Run an estep (variant calling step) job per section in the splitList
process CAVEMAN_ESTEP {
    script:
    """
    ...
    -process estep \
    –index ${splitlist_channel; 1,2,3...n.}
    """
}

// Merge the CaVEMan results into a single VCF result file
process CAVEMAN_RESULTS {
    script:
    """
    ...
    -process merge_results \
    –index 1
    """
}

// Add UUIDs to each VCF entry
process CAVEMAN_ID {
    script:
    """
    ...
    -process add_ids \
    –index 1
    """
}

// Flag (post process) the VCF results
process CAVEMAN_FLAG {
    script:
    """
    ...
    -process flag \
    –index 1
    """
}