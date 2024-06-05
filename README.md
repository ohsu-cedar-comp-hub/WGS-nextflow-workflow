# WGS Nextflow Workflow

This Nextflow workflow aligns, performs somatic variant calls, and annotates variants on whole genome sequencing matched tumor-normal data. The pipeline is intended to run with a container (Docker or Singularity). Submodules of the workflow are described below. A parameter file is passed as input using `-params-file <my-params.json>`.

## Getting Started

### Running with Docker/Singularity

Requires install of [Singularity](https://sylabs.io/docs/) or [Docker](https://www.docker.com/get-started/).

Check that Java version is 11 through 21. Nextflow has been checked to run with Java 11-21. Otherwise install Java 17 via SDKMAN shown [here](https://www.nextflow.io/docs/latest/install.html).  

```
java -version
```

Install [Nextflow](https://www.nextflow.io/docs/latest/install.html).

```
# Install Nextflow
curl -s https://get.nextflow.io | bash

# Make Nextflow executable
chmod +x nextflow

# Move Nextflow into an executable path:
sudo mv nextflow /usr/local/bin

# Confirm install:
nextflow info

# Update to latest Nextflow version
nextflow self-update
```

Dockerfiles for building images are provided for each tool in the pipeline in the tool directory, `tools/<tool name>/Dockerfile`:
- fastqc 
- multiqc
- trimmomatic
- bwa
- samtools
- gatk
- bcftools
- snpeff

Specify paths to your containers in the [params file](example/nf_channels_params_template.json).  

Singularity images can be pulled from Docker images. [Read more about Singularity here.](https://docs.sylabs.io/guides/3.2/user-guide/cli/singularity_pull.html)

Enable running with your preferred container in your Nextflow [config file](https://www.nextflow.io/docs/latest/config.html).  
Here is an example with Singularity. You will need to enable auto mounts. 

```
singularity {
    enabled = true
    autoMounts = true
    cacheDir = '<your cache dir>'
        }
```

Invoking the workflow   

When running Nextflow, you will need:
- the script `(workflows/<step>/<script name>.nf)`
- the parameters file (an example of the params file is available in this repo [here](example/nf_channels_params_template.json).
- the nextflow config file


```
nextflow run <nextflow script>.nf \
-params-file <parameters>.json \
-c nextflow.config
```

### Data setup

Example Ceph Bucket Config

```
aws {
    accessKey='<KEY>'
    secretKey='<SECRET>'
    client {
        endpoint = 'https://rgw.ohsu.edu'
        s3PathStyleAccess = true
    }
}
```


## Workflow design

**Pre-alignment QC**

1. Initial QC with FastQC   
2. Trim with Trimmomatic  
3. fastQC on trimmed reads  
4. MultiQC on all fastQC output files  

**Alignment**  

5. Alignment with bwa-mem2  
6. Sort with samtools  
7. Mark duplicates with gatk MarkDuplicates  
8. Sort and index with samtools

**Somatic Variant Calling**  

9. Get pileup summaries with gatk GetPileupSummaries
10. Calculate contamination with gatk CalculateContamination
11. Somatic variant calling with gatk Mutect2, parallelized per chromosome
12. Process Mutect2 output  
    12a. bgzip each of the file outputs  
    12b. Concatenate into single VCF with bcftools concat  
    12c. Normalize VCF with bcftools normalize  
    12d. Merge stats files  
    12e. Merge f1r2 files  
13. Learn technical artifact prior probability with gatk LearnReadOrientationModel
14. Filter variant calls with gatk Mutect2

**Annotation**  

15. Annotate variants with snpEff   


## Alignment workflow

Description of the [WGS nextflow alignment pipeline](https://github.com/ohsu-cedar-comp-hub/WGS-nextflow-workflow) with a comparison to the [GDC pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/).

![align_pipeline](https://github.com/ohsu-cedar-comp-hub/WGS-nextflow-workflow/assets/136844363/038d08a6-1621-43e1-b8f9-da62a89871f8)

### Tools 

- Nextflow (DSL2)
- Singularity
- fastqc v0.11.7
- multiqc v1.21
- trimmomatic 0.38
- bwa-mem2 v2.2.1
- samtools 1.19
- gatk4 4.4.0.0
  - Picard MarkDuplicates

|Step | Parameter | COH WGS nextflow pipeline | GDC pipeline |
| --- | --------- | ------------------------- | ------------ |
| 0 | Reference genome | GRCh38.d1.vd1 | GRCh38.d1.vd1 |
| 1 | BWA algorithm | bwa-mem2 | bwa-mem |
| 2 | BAM sort | samtools | picard SortSam |
| 3 | Mark duplicates | picard MarkDuplicates | picard MarkDuplicates |
| 4 | Co-cleaning workflow | NA | gatk RealignerTargetCreator, IndelRealigner, BaseRecalibrator, PrintReads | 

## Usage  

Command line call:
```
nextflow run workflows/alignment/align_wf.nf -c nextflow.config -params-file <your_params_file>
```

### 0. Pre-alignment: quality filtering and trimming  

**Run fastQC**  
**File input**: Raw paired end reads in fastq file format  
**File output**: .html and .zip fastQC analysis files
fastQC requires a directory to be made before running. Make a directory in the same output path as specified in the parameters file.  

**Trimmomatic**  
After QC, trimmomatic is run in paired-end mode.   
**File input**: Raw paired end reads in fastq format and a file of sequencing adapters used in the library prep (TruSeq3-PE.fa)  
**File output**: Two fastq files of quality filtered, adapter-trimmed paired end reads (trim read 1 and trim read 2), and two fastq files of unpaired reads.

**Script in nextflow:**
```Shell
java -jar /bin/trimmomatic.jar \
    PE -phred33 \
    ${reads[0]} \
    ${reads[1]} \
    ${reads[0].simpleName}_1P.fastq.gz \
    ${reads[0].simpleName}_1U.fastq.gz \
    ${reads[1].simpleName}_2P.fastq.gz \
    ${reads[1].simpleName}_2U.fastq.gz \
    ILLUMINACLIP:"${truseq3pefile}":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
```

*Note: Currently, adjusting arguments for trimmomatic is not supported with parameters and requires manually editing the script block, but is noted for future development*

**fastQC on trimmed reads**  
After running trimmomatic, run fastQC again on trimmed reads.  
**File input**: Trimmed, quality filtered fastq files of paired end reads (output from trimmomatic)  
**File output**: .html and .zip files of fastQC analysis

**MultiQC**  
**File input**: fastQC results from raw reads and trimmed reads    
**File output**: .html multiQC report

### 1. Alignment with BWA-MEM2
 
**File input:** Trimmed, quality filtered, paired-end .fastq files, and an indexed reference genome fasta file.   
**File output:** Aligned, unsorted BAM file   
 
**Script in nextflow:**
``` Shell
bwa-mem2 mem \
    -K 100000000 -t 6 -Y -M \
    -R "@RG\\tID:${params.ID}\\tLB:no_library\\tPL:illumina\\tPU:none\\tSM:${trim_read1.simpleName}" \
    ${params.idx} ${trim_read1} ${trim_read2} |
    samtools view -Sb -@ 4 > ${trim_read1.simpleName}.bam
```
**Comparison GDC script**:  
```Shell
bwa mem \
    -t 8 \
    -T 0 \
    -R <read_group> \
    <reference> \
    <fastq_1.fq.gz> \
    <fastq_2.fq.gz> |
    samtools view \
    -Shb -o <output.bam>
```

### 2. Sort and index
**File input**: Unsorted BAM file  
**File output**: Sorted and indexed BAM file (.bam file and .bam.bai file)  

**Script in nextflow:**
```Shell
samtools sort ${bam_unsorted} > ${bam_unsorted.baseName}_sorted_indexed.bam
samtools index ${bam_unsorted.baseName}_sorted_indexed.bam > ${bam_unsorted.baseName}_sorted_indexed.bam.bai
```  

**Comparison GDC script**:   
```Shell
java -jar picard.jar SortSam \
    CREATE_INDEX=true \
    INPUT=<input.bam> \
    OUTPUT=<output.bam> \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=STRICT
```

### 3. Mark duplicates
**File input:** Sorted BAM file   
**File output:** Unsorted BAM file of marked duplicates and a metrics txt file  

**Script in nextflow:**  
```Shell
gatk MarkDuplicates \
    -I ${bam_sorted} \
    -O ${bam_sorted.baseName}_marked_duplicates.bam \
    -M ${bam_sorted.baseName}_marked_duplicates_metrics.txt \
    --VALIDATION_STRINGENCY LENIENT
```  

**Comparison GDC script**:  
```Shell
java -jar picard.jar MarkDuplicates \
    CREATE_INDEX=true \
    INPUT=<input.bam> \
    VALIDATION_STRINGENCY=STRICT
```
 
### 4. Sort 
**File input**: Unsorted duplicate-marked BAM file  
**File output**: Sorted duplicate-marked BAM file

Script in Nextflow:
```Shell
samtools sort ${bam_duplicates_unsorted} > ${bam_duplicates_sorted.baseName}_sorted.bam
```
Comparison GDC script:
```Shell
java -jar picard.jar SortSam \
    I=input.bam \
    O=sorted.bam \
    SORT_ORDER=coordinate
```

## Variant Calling Workflow and Annotation Workflow


### Tools
- Nextflow (DSL2)
- Singularity or Docker
- htslib/1.2
  - bgzip
- bcftools/1.12
  - concat
  - normalize
  - merge
- GATK:4.4.0.0
  - GeneratePileUpSummaries
  - CalculateContamination
  - Mutect2
  - LearnReadOrientation 
  - MergeMutectStats
  - FilterMutectCalls
- SNPEff/5.2c


## Usage

Command line call:
```
nextflow run workflows/calling/svc_wf.nf -c nextflow.config -params-file <your_params_file>
```

### 1. GetPileupSummaries
*Generate read counts for estimating of contamination in the tumor tissue using GATK.*  
**File input**: Duplicate-marked BAM files of tumor and matched normal and common germline variants VCF (exac file)  
**File output**: Pileup summary table of ref and alt counts, and allele frequencies at each position

### 2. Calculate Contamination
*Calculate fraction of normal cell contaminants in tumor sample*  
**File input**: Tumor and normal pileup summary tables  
**File output**: Contamination table

### 3. Run Mutect2 
*Call somatic variants*  
**File input**: Duplicate-marked BAM files of tumor and matched normal, reference genome, germline resource VCF (exac), panel of normals (PoN) VCF  
**File output**: 23 f1r2.tar.gz files (1 per chromosomes 1-22 and X), 23 unfiltered tumor VCFs, and 23 VCF stats files

Script in nextflow:  
```
gatk Mutect2 \
    -R ${params.mutect_idx} \
    -I ${tumor_bam_sorted} \
    -I ${normal_bam_sorted} \
    --panel-of-normals ${params.pon} \
    -L ${chrom} \
    --germline-resource ${params.gnomad} \
    -O ${sample_id}_${chrom}_unfiltered.vcf \
    --f1r2-tar-gz ${sample_id}_${chrom}_f1r2.tar.gz
```

### 4. Organize Mutect2 output 
#### A. bgzip each of the 23 files
#### B. concatenate 23 chromosome-specific bgzipped files into one VCF file  
#### C. normalize VCF file  
#### D. merge 23 chromosome-specific stats files into one stats file  
#### E. combine 23 chromosome-specific f1r2.tar.gz files into one file  

### 5. GATK LearnReadOrientationModel
_Learn the read orientation model to refine variant calls by removing technical artifacts._   
**File input**: 23 individual f1r2 files (one per chromosome), containing the outputs from Mutect2  
**File output**: artifacts priors.tar.gz file    

### 6. GATK FilterMutectCalls
_Apply filters to Mutect2 variant calls_  
**File input**: Unfiltered, combined across chromosomes normalized VCF and VCF stats file, genome reference (including index & dict), f1r2 read orientation model and contamination table.  
**File output**: Filtered VCF   

Script in nextflow:   

```
gatk FilterMutectCalls \
    -R ${params.mutect_idx} \
    -V ${unfiltered_vcf} \
    --tumor-segmentation ${segmentation_table} \
    --contamination-table ${contamination_table} \
    --ob-priors ${read_orientation_model} \
    -O ${sample_id}_filtered.vcf \
    --stats ${vcf_stats}
```

### 7. SNPEff Annotation
_Annotate variants using SNPEff_  
**File input**: Filtered VCF and genome version GRCh38.86 from snpEff pre-built database  
**File output**: Annotated variants VCF

### Output directory structure 

Running the pipeline generates a total of 97 files in 14 directories, structured like this:
```
├── aligned
│   ├── YOURSAMPLEID_G.bam
│   ├── YOURSAMPLEID_T.bam
│   ├── markduplicates
│   │   ├── YOURSAMPLEID_G_sorted_marked_duplicates.bam
│   │   ├── YOURSAMPLEID_G_sorted_marked_duplicates_metrics.txt
│   │   ├── YOURSAMPLEID_T_sorted_marked_duplicates.bam
│   │   ├── YOURSAMPLEID_T_sorted_marked_duplicates_metrics.txt
│   │   └── sorted
│   │       ├── YOURSAMPLEID_G_sorted_marked_duplicates_sorted.bam
│   │       ├── YOURSAMPLEID_G_sorted_marked_duplicates_sorted.bam.bai
│   │       ├── YOURSAMPLEID_T_sorted_marked_duplicates_sorted.bam
│   │       └── YOURSAMPLEID_T_sorted_marked_duplicates_sorted.bam.bai
│   └── sort_index
│       ├── YOURSAMPLEID_G_sorted.bam
│       └── YOURSAMPLEID_T_sorted.bam
├── fastqc
│   ├── YOURSAMPLEID_G_R1_1P_fastqc.html
│   ├── YOURSAMPLEID_G_R1_1P_fastqc.zip
│   ├── YOURSAMPLEID_G_R1_fastqc.html
│   ├── YOURSAMPLEID_G_R1_fastqc.zip
│   ├── YOURSAMPLEID_G_R2_2P_fastqc.html
│   ├── YOURSAMPLEID_G_R2_2P_fastqc.zip
│   ├── YOURSAMPLEID_G_R2_fastqc.html
│   ├── YOURSAMPLEID_G_R2_fastqc.zip
│   ├── YOURSAMPLEID_T_R1_1P_fastqc.html
│   ├── YOURSAMPLEID_T_R1_1P_fastqc.zip
│   ├── YOURSAMPLEID_T_R1_fastqc.html
│   ├── YOURSAMPLEID_T_R1_fastqc.zip
│   ├── YOURSAMPLEID_T_R2_2P_fastqc.html
│   ├── YOURSAMPLEID_T_R2_2P_fastqc.zip
│   ├── YOURSAMPLEID_T_R2_fastqc.html
│   └── YOURSAMPLEID_T_R2_fastqc.zip
├── filtered
│   └── YOURSAMPLEID_filtered.vcf
├── multiqc
│   └── multiqc_report.html
├── summaries
│   ├── YOURSAMPLEID_G_sorted_marked_duplicates_sorted.getpileupsummaries.table
│   └── YOURSAMPLEID_T_sorted_marked_duplicates_sorted.getpileupsummaries.table
├── svc
│   ├── annotated_variants
│   │   └── YOURSAMPLEID_annotated_variants.vcf
│   ├── YOURSAMPLEID_chr10_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr10_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr11_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr11_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr12_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr12_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr13_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr13_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr14_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr14_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr15_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr15_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr16_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr16_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr17_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr17_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr18_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr18_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr19_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr19_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr1_unfiltered.vcf
│   ├── YOURSAMPLEID_chr1_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr1_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr20_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr20_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr21_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr21_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr22_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr22_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr2_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr2_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr3_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr3_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr4_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr4_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr5_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr5_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr6_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr6_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr7_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr7_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr8_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr8_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chr9_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chr9_unfiltered.vcf.gz.tbi
│   ├── YOURSAMPLEID_chrX_unfiltered.vcf.gz
│   ├── YOURSAMPLEID_chrX_unfiltered.vcf.gz.tbi
│   ├── f1r2files
│   │   └── YOURSAMPLEID_read-orientation-model.tar.gz
│   └── sort_index
│       ├── YOURSAMPLEID_normalized.vcf.gz
│       ├── YOURSAMPLEID_normalized.vcf.gz.tbi
│       ├── YOURSAMPLEID_unfiltered_sorted.vcf.gz
│       ├── YOURSAMPLEID_unfiltered_sorted.vcf.gz.tbi
│       ├── YOURSAMPLEID_unfiltered.vcf.all.stats
│       └── YOURSAMPLEID_unfiltered.vcf.gz
├── tables
│   ├── YOURSAMPLEID_T_sorted_marked_duplicates_sorted_contamination_table
│   └── YOURSAMPLEID_T_sorted_marked_duplicates_sorted_segmentation_table
└── trim_reads
    ├── YOURSAMPLEID_G_R1_1P.fastq.gz
    ├── YOURSAMPLEID_G_R1_1U.fastq.gz
    ├── YOURSAMPLEID_G_R2_2P.fastq.gz
    ├── YOURSAMPLEID_G_R2_2U.fastq.gz
    ├── YOURSAMPLEID_T_R1_1P.fastq.gz
    ├── YOURSAMPLEID_T_R1_1U.fastq.gz
    ├── YOURSAMPLEID_T_R2_2P.fastq.gz
    └── YOURSAMPLEID_T_R2_2U.fastq.gz
```

*Note: future development will reduce the number of files in intermediary steps that are automatically output by the pipeline*