Current workflow design for single and paired FASTQ files: 
  1.) QC with Fastqc
  2.) QC with FastQC and MultiQC
  2.) Trim with Trimmomatic
  3.) Align with BWA-MEM2

Prerequisits: environment set-up
  1.) Conda create -n NAME python=3.11
  2.) Install nextflow (conda install nextflow)
  2.) Install fastqc (conda install fastqc)
  3.) Install multiqc (pip install multiqc)
  4.) Install trimmomatic (conda install trimmomatic) 
  5.) Install bwa-mem2 (conda install bwa-mem2) 
  6.) Install sam-tools (conda install -c bioconda samtools)

