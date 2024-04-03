Current workflow design for single and paired FASTQ files: 
  1.) QC with Fastqc
  2.) Trim with TrimmomaticPE
  3.) Fastqc on trimmed reads
  4.) MultiQC on all fastqc output files
  5.) Alignment with bwa-mem2
  6.) sort and index with samtools 
  7.) Mark duplicates with gatk MarkDuplicates
  8.) somatic variant calling with gatk Mutect2
  9.) Annotating variants with snpEff

Prerequisites: environment set-up
  1.) Conda create -n NAME python=3.11
  2.) Install nextflow (conda install nextflow)
  2.) Install fastqc (conda install fastqc)
  3.) Install multiqc (pip install multiqc)
  4.) Install trimmomatic (conda install trimmomatic) 
  5.) Install bwa-mem2 (conda install bwa-mem2) 
  6.) Install sam-tools (conda install -c bioconda samtools)
  7.) Install gatk (conda install gatk) or download from gatk github
  8.) Install snpEff (conda install snpEff) or download from snpEff source page 
      - download reference genome with snpEff download GRCh38.86 

