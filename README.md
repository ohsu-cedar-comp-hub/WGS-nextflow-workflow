# WGS Nextflow Workflow

This Nextflow workflow aligns, performs somatic variant calls, and annotates variants on whole genome sequencing matched tumor-normal data. The pipeline can be run with conda or Docker/Singularity, but it is recommended to use the containerized approach with Singularity. Submodules of the workflow are described below and can be run independently assuming the necessary input files exist. A parameter file is passed as input using -params-file <my-params.json>, which can be generated using the templating script.

## Getting Started

### Running with Anaconda
Requires install of [Anaconda](https://docs.anaconda.com/free/anaconda/install/index.html)

Create conda environment
```
conda create -n wgs python=3.11
```

Activate conda environment
```
conda activate wgs
```

Install dependencies
```
conda install -c bioconda nextflow fastqc trimmomatic bwa-mem2 samtools gatk snpEff --yes
pip install multiqc
```

Download reference genome file with snpEff download GRCh38.86

Check that Java version is 11 through 21. Nextflow has been checked to run with Java 11-21. Otherwise install Java 17 via SDKMAN shown [here](https://www.nextflow.io/docs/latest/install.html)  

```
java -version
```

Test Nexflow is successfully installed
```
nextflow info
```

Invoke workflow  

When running Nextflow, you will need:  
- the script `(workflows/<step>/<script name>.nf)`  
- the parameters file (an example of the params file and a python script to generate it is available in this repo here: `/config/example_pair_params.json`) 
- the nextflow [config file](https://www.nextflow.io/docs/latest/config.html)  
  
```
nextflow run <nextflow script>.nf \
-params-file <parameters>.json \
-c nextflow.config
```

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
- mutect *[RL] EDIT: redundant with gatk?
- bcftools
- snpeff

Add Singularity to your Nextflow config file

```
singularity {
    enabled = true
    autoMounts = true
    cacheDir = '<your cache dir>'
        }
```
Invoke Workflow   

When running Nextflow, you will need:
- the script `(workflows/<step>/<script name>.nf)`
- the parameters file (an example of the params file and a python script to generate it is available in this repo here: `/config/example_pair_params.json`)
- the nextflow [config file](https://www.nextflow.io/docs/latest/config.html)
- the path to your singularity image (.sif file)  
<br>

```
nextflow run <nextflow script>.nf \
-params-file <parameters>.json \
-c nextflow.config \
-with-singularity /path/to/<name>.sif
```


---
[RL] EDIT **Move exacloud-specific steps to github submodule**

    ### Running on Exacloud [internal]

    This workflow was built to run on the [Exacloud HPC](https://wiki.ohsu.edu/display/ACC/Exacloud) and installations or module loading for the prerequisites are specific to this environment.

    **Nextflow setup**

    ```
    ## Create conda environment and install python
    conda create -n nextflow_env python=3.11
    ## activate environment
    conda activate nextflow_env
    ## install nextflow
    conda install nextflow 
    ## self update nextflow to newest release
    nextflow self-update
    ```

    **Singularity setup**  
    See more on running Singularity on Exacloud [here](https://wiki.ohsu.edu/display/ACC/Exacloud%3A+Singularity).

    ```
    ## Load singularity module 
    module load /etc/modulefiles/singularity/current

    ## singularity cache directory on gscratch
    export SINGULARITY_CACHEDIR=/home/exacloud/gscratch/<group>/<username>/singularity

    ## remote login to quay
    singularity remote login --username myuser docker://quay.io

    ## set temporary directory
    export TMPDIR=`pwd`/tmp
    ```

    Pull Singularity images

    ```
    cd <directory where your .sif files will live>
    singularity pull <name>.sif docker://quay.io/ohsu-comp-bio/<name>
    ```
    
    There are currently separate images for each tool. You will need to build the following containers using the Dockerfile associated with each tool and convert them to Singularity image files if working on a HPCC:
    - FastQC 
    - MultiQC
    - trimmomatic
    - bwa
    - samtools
    - GATK4
    - Mutect2
    - bcftools
    - snpeff


    Example:

    ```
    cd config_files
    singularity pull bwa.sif docker://quay.io/ohsu-comp-bio/bwa
    ```

    **environment summary** 
    Nextflow version 23.10.1  
    Singularity version 3.8.0-1.el7

    **End of section to put into submodule**
---

### Data setup (rename?)

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
6. Sort and index with samtools  
7. Mark duplicates with gatk MarkDuplicates  
8. Sort and index with samtools * EDIT [RL]: is it possible to condense these 2 sorting and indexing steps?  

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


## Alignment workflow [WIP]

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

### 0. Pre-alignment: quality filtering and trimming  

**Run fastQC**  
**File input**: Raw paired end reads in fastq file format  
**File output**: .html and .zip fastQC analysis files
fastQC requires a directory to be made before running. Make a directory in the same output path as specified in the parameters file.  


```
## make output directory for fastqc
mkdir output/fastqc

## run nextflow
nextflow run workflows/qc/fastqc.nf \
-params-file <params-file>.json \
-c nextflow.config \
-with-singularity fastqc.sif 
```

**Trimmomatic**  
After QC, trimmomatic is run in paired-end mode.   
**File input**: Raw paired end reads in fastq format and a file of sequencing adapters used in the library prep (TruSeq3-PE.fa)  
**File output**: Two fastq files of quality filtered, adapter-trimmed paired end reads (trim read 1 and trim read 2), and two fastq files of unpaired reads.

```
nextflow run trimmomatic.nf \
—params-file <params-file>.json \
-c nextflow.config \
-with-singularity trimmomatic.sif
```
**Script in nextflow:**
```Shell
java -jar /bin/trimmomatic.jar \
    PE -phred33 \
    ${read1} \
    ${read2} \
    ${read1.baseName}_1P.fastq.gz \
    ${read1.baseName}_1U.fastq.gz \
    ${read2.baseName}_2P.fastq.gz \
    ${read2.baseName}_2U.fastq.gz \
    ILLUMINACLIP:"${truseq3pefile}":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
```
  
**fastQC on trimmed reads**  
After running trimmomatic, run fastQC again on trimmed reads.  
**File input**: Trimmed, quality filtered fastq files of paired end reads (output from trimmomatic)  
**File output**: .html and .zip files of fastQC analysis

```
nextflow run workflows/qc/trim_fastqc.nf \
-params-file <params-file>.json \
-c nextflow.config \
-with-singularity fastqc.sif 
```

**MultiQC**  
**File input**: fastQC results from raw reads and trimmed reads    
**File output**: .html multiQC report
```
nextflow run workflows/qc/multiqc.nf \
-params-file <params-file>.json \
-c nextflow.config \
-with-singularity multiqc.sif 
```

### 1. Alignment with BWA-MEM2
 
**File input:** Trimmed, quality filtered, paired-end .fastq files, and an indexed reference genome fasta file.   
**File output:** Aligned, unsorted BAM file   
 
```
nextflow run bwamem2.nf \
—params-file <params-file>.json \
-c nextflow.config \
-with-singularity bwa.sif>
```

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
    -Shb -o <output.bam> -
```

### 2. Sort and index
**File input**: Unsorted BAM file  
**File output**: Sorted and indexed BAM file (.bam file and .bam.bai file)  

```
nextflow run sort_and_index.nf \
—params-file <params-file>.json \
-c nextflow.config \
-with-singularity samtools.sif>
```

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
 
```
nextflow run mark_duplicates.nf \
—params-file <params-file>.json \
-c nextflow.config \
-with-singularity gatk.sif
```

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

```
nextflow run sort_marked_duplicates.nf \
—params-file <params-file>.json \
-c nextflow.config \
-with-singularity samtools.sif
```

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

## Variant Calling Workflow and Annotation Workflow [WIP] 


![variantcall_pipeline](https://github.com/ohsu-cedar-comp-hub/WGS-nextflow-workflow/assets/136844363/1aa03724-3ba7-466d-ab0e-de90ca1d56a3)

### Tools
- Nextflow (DSL2)
- Singularity 
- bcftools/1.10.2
  - concat
  - normalize
  - merge
- GATK:4.4.0.0
  - GeneratePileUpSummaries
  - CalculateContamination
  - LearnReadOrientation 
  - Mutect2
- SNPEff/5.2c
- htslib/1.10.2

## Usage

### 1. GetPileupSummaries
*Generate read counts for estimating of contamination in the tumor tissue using GATK.*  
**File input**: Duplicate-marked BAM files of tumor and matched normal and common germline variants VCF (exac file)  
**File output**: Pileup summary table of ref and alt counts, and allele frequencies at each position

```
nextflow run get_pileup_summaries.nf \
—params-file <params-file>.json \
-c nextflow.config \
-with-singularity gatk.sif
```

### 2. Calculate Contamination
*Calculate fraction of normal cell contaminants in tumor sample*  
**File input**: Tumor and normal pileup summary tables  
**File output**: Contamination table

```
nextflow run calculate_contamination.nf \
—params-file <params-file>.json \
-c nextflow.config \
-with-singularity gatk.sif
```

### 3. Run Mutect2 
*Call somatic variants*  
**File input**: Duplicate-marked BAM files of tumor and matched normal, reference genome, germline resource VCF (exac), panel of normals (PoN) VCF  
**File output**: 23 f1r2.tar.gz files (1 per chromosomes 1 -22 and X), 23 unfiltered tumor VCFs, and 23 VCF stats files

```
nextflow run mutect2.nf \
—params-file <params-file>.json \
-c nextflow.config \
-with-singularity mutect.sif
```
Script in nextflow:  
```
gatk Mutect2 \\
    -R ${mutect_idx} \\
    -I ${tumor_bam_sorted} \\
    -I ${normal_bam_sorted} \\
    -normal ${normal_bam_sorted.baseName} \\
    --panel-of-normals ${pon} \\
    --germline-resource ${germline_resource} \\
    -L ${chrom} \\ filesChannel() containing array of chromsomes, chr1 ... chrX 
    -O ${tumor_bam_sorted.baseName}_${chrom}_unfiltered.vcf \\
    --f1r2-tar-gz ${tumor_bam_sorted.baseName}_${chrom}_f1r2.tar.gz \\
    -stats ${tumor_bam.baseName}_${chrom}_unfiltered.vcf.stats
```

```docker run quay.io://ohsu-comp-bio/gatk:4.4.0.0 gatk Mutect2 -R reference.fa -I tumor.bam -I normal.bam -normal normal_name --intervals chrom -pon gnomad_panel_of_normals.vcf -germline-resource af-only_exac.vcf --f1r2-tar-gz "${file%.bam}.f1r2.tar.gz -O "${file%.bam}.unfiltered.vcf"```

### 4. Process Mutect2 output 
#### A. bgzip each of the 23 files 

#### B. concatenate 23 chromosome-specific bgzipped files into one VCF file
#### C. normalize VCF file
#### D. merge 23 chromosome-specific stats files into one stats file
#### E. combine 23 chromosome-specific f1r2.tar.gz files into one file

### 5. GATK LearnReadOrientationModel
_Learn the read orientation model to refine variant calls by removing technical artifacts._ 
**File input**: f1r2 file, containing the combined (per-chromosome) outputs from Mutect2
**File output**: artifacts-priors.tar.gz file  

### 6. GATK FilterMutectCalls
_Apply filters to Mutect2 variant calls_  
**File input**: Unfiltered, combined across chromosomes normalized VCF and VCF stats file, genome reference (including index & dict), f1r2 read orientation model and contamination table.  
**File output**: Filtered VCF   

```
nextflow run annotate_variants.nf \
—params-file <params-file>.json \
-c nextflow.config \
-with-singularity mutect.sif
```
Script in nextflow:   

```
    gatk FilterMutectCalls \
    -R ${mutect_idx} \
    -V ${unfiltered_vcf} \
    --contamination-table ${contamination_table} \
    -O ${unfiltered_vcf.baseName}_filtered.vcf \
    --read-index ${mutect_idx_fai} \
    --sequence-dictionary ${mutect_dict} \
    --ob-priors ${prior-artifacts.tar.gz} \
    --stats ${vcf_stats}
```

```docker run quay.io://ohsu-comp-bio/gatk:4.4.0.0 gatk FilterMutectCalls -R reference.fa -V  "${file}.unfiltered.normalized.vcf" -contamination-table contamination.table --orientation-bias-artifact-priors read-orientation-model.tar.gz -O "${file%.unfiltered.normalized.vcf.gz}.filtered.vcf" ```

### 7. SNPEff Annotation
_Annotate variants using SNPEff_
**File input**: Filtered VCF and genome version GRCh38.86 from snpEff pre-built database
**File output**: Annotated variants VCF

```
nextflow run annotate_variants.nf \
—params-file <params-file>.json \
-c nextflow.config \
-with-singularity snpeff.sif
```


