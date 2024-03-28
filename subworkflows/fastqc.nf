// Define the process for running FastQC
process fastQC {
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    path read1
    path read2
    
    script:
    """
     fastqc -o ${params.outdir}/fastqc ${read1}  
     fastqc -o ${params.outdir}/fastqc ${read2} 
    """
}
