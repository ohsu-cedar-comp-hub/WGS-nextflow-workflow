# WGS Nextflow Workflow

This Nextflow workflow aligns, calls, and annotates whole genome sequencing data using Singularity. Submodules of the workflow are described and can be run independently.

Below is the current alignment and variant calling workflow design for matched tumor-normal data generated by whole genome sequencing.  

## Prerequisites:

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
  
There are currently separate images for each tool. You will need to pull the following from quay.io/ohsu-comp-bio:
- fastqc [wip]
- trimmomatic
- bwa
- samtools
- gatk
- mutect
- bcftools
- snpeff [wip]


Example:

```
cd config_files
singularity pull bwa.sif docker://quay.io/ohsu-comp-bio/bwa
```

Add Singularity to your Nextflow config file

```
singularity {
    enabled = true
    autoMounts = true
    cacheDir = '<your cache dir>'
    container = ‘file://path/to/<name>.sif'
        }
```

Add environmental paths to your config file

for fastqc:
```
env {
    PATH = "$PATH:/usr/local/FastQC/"
}
```

Run Nextflow with Singularity  
When running Nextflow, you will need:
- the script `(workflows/<step>/<script name>.nf)`
- the parameter file (an example of the parameter file and a python script that can be used to generate one's own unique parameter file are available in this repo at `WGS-nextflow-workflow/example/normal_params_file_template.json` and `WGS-nextflow-workflow/example/render_normal_params.py`)
- the nextflow [config file](https://www.nextflow.io/docs/latest/config.html)
- the path to your singularity image (.sif file)  
<br>

```
nextflow run <nextflow script>.nf \
-params-file <parameters>.json \
-c nextflow.config \
-with-singularity /path/to/<name>.sif
```

**Environment Summary**

Nextflow version 23.10.1  
Singularity version 3.8.0-1.el7

## Quick run 

Activate nextflow environment and load singularity

**Pre-alignment QC**   
1.) QC with fastQC 

 ```
## make output directory for fastqc
mkdir output/fastqc
## run nextflow
nextflow run workflows/qc/fastqc.nf -params-file <params-file>.json -c nextflow.config -with-singularity fastqc.sif 
``` 
  
add file outputs from fastQC to params file  
  
2.) Trim with TrimmomaticPE  
`nextflow run trimmomatic.nf -params-file <params-file>.json -c nextflow.config -with-singularity trimmomatic.sif`  
  
add file outputs from trimmomatic to params file  
  
3.) fastQC on trimmed reads   
`nextflow run workflows/qc/trim_fastqc.nf -params-file <params-file>.json -c nextflow.config -with-singularity fastqc.sif`  
  
add file outputs from fastqc to params file  
  
4.) MultiQC on all fastQC output files   
`nextflow run workflows/qc/multiqc.nf -params-file <params-file>.json -c nextflow.config -with-singularity fastqc.sif`  

**Alignment**  
5.) Alignment with bwa-mem2  
`nextflow run workflows/bwa/bwamem2.nf -params-file <params-file>.json -c nextflow.config -with-singularity bwa.sif`

add unsorted bam output to params file

6.) Sort and index with samtools  
`nextflow run workflows/samtools/sort_and_index.nf -params-file <params-file>.json -c nextflow.config -with-singularity samtools.sif`  

add sorted indexed bam output to params file  

7.) Mark duplicates with gatk MarkDuplicates  
`nextflow run workflows/gatk/mark_duplicates.nf -params-file <params-file>.json -c nextflow.config -with-singularity gatk.sif`  
  
**Somatic Variant Calling**  
8.) Somatic variant calling with gatk Mutect2   
`nextflow run workflows/mutect/mutect2.nf -params-file <params-file>.json -c nextflow.config -with-singularity mutect.sif`  
<< EDIT not sure why we have separate mutect and gatk images?

add f1r2 file to params file  
add unfiltered vcf (svc vcf output from mutect) to params file

9.) Get segmentation and contamination tables  
`nextflow run workflows/gatk/calculate_contamination.nf -params-file <params-file>.json -c nextflow.config -with-singularity gatk.sif`  

*split by chromosome with bash script here* 

10.) Learn orientation bias model  
`nextflow run workflows/gatk/learn_read_orientation_model.nf -params-file <params-file>.json -c nextflow.config -with-singularity gatk.sif`
add output to params file

11.) Get pileup summaries  
`nextflow run workflows/gatk/get_pileup_summaries.nf -params-file <params-file>.json -c nextflow.config -with-singularity gatk.sif`  
add output file to params file  

12.) Sort, index, normalize with bcftools  
`nextflow run workflows/bcftools/vcf_sort_index_normalize.nf -params-file <params-file>.json -c nextflow.config -with-singularity bcftools.sif`  
add output to params file

13.) Filter variant calls with gatk Mutect2  
`nextflow run workflows/mutect/filter_mutect.nf -params-file <params-file>.json -c nextflow.config -with-singularity mutect.sif`  
add filtered vcf to params file 

**Annotation**

14.) Annotate variants with snpEff   
`nextflow run workflows/snpeff/annotate_variants.nf -params-file <params-file>.json -c nextflow.config -with-singularity snpeff.sif`  

15.) Calculate allele frequencies (?)

## Workflow description

### Alignment workflow [WIP]

Description of the [WGS nextflow alignment pipeline](https://github.com/ohsu-cedar-comp-hub/WGS-nextflow-workflow) with a comparison to the [GDC pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/).

![align_pipeline](https://github.com/ohsu-cedar-comp-hub/WGS-nextflow-workflow/assets/136844363/038d08a6-1621-43e1-b8f9-da62a89871f8)

### Tools 

- Nextflow (DSL2)
- Singularity
- trimmomatic 0.38
- bwa-mem2 {VERSION} << we are cloning the bwa-mem2 repo for this, for some reason I can't find the command to run to check version
- samtools 1.19
- gatk4 4.4.0.0

|Step | Parameter | COH WGS nextflow pipeline | GDC pipeline |
| --- | --------- | ------------------------- | ------------ |
| 0 | Reference genome | GRCh38.d1.vd1 | GRCh38.d1.vd1 |
| 0 | Read groups | -R flag | Aligned separately then merged << Unsure about this section |
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

Ensure your path is updated in the nextflow config file. 

```
env {
    PATH = "$PATH:/usr/local/FastQC/"
}
```
  
```
## make output directory for fastqc
mkdir output/fastqc
## run nextflow
nextflow run workflows/qc/fastqc.nf -params-file <params-file>.json -c nextflow.config -with-singularity fastqc.sif 
```

**Trimmomatic**  
After QC, trimmomatic is run in paired-end mode.   
**File input**: Raw paired end reads in fastq format and a file of sequencing adapters used in the library prep (TruSeq3-PE.fa)  
**File output**: Two fastq files of quality filtered, adapter-trimmed paired end reads (trim read 1 and trim read 2), and two fastq files of unpaired reads.

```
nextflow run trimmomatic.nf —params-file <my-params.json> -c <my-nextflow.config> -with-singularity <image.sif>
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
nextflow run workflows/qc/trim_fastqc.nf -params-file <params-file>.json -c nextflow.config -with-singularity fastqc.sif 
```

**MultiQC**  
**File input**: fastQC results from raw reads and trimmed reads    
**File output**: .html multiQC report
```
nextflow run workflows/qc/multiqc.nf -params-file <params-file>.json -c nextflow.config -with-singularity fastqc.sif 
```

### 1. Alignment with Burrows-Wheeler Aligner (BWA)
 
**File input:** Trimmed, quality filtered, paired-end .fastq files, and an indexed reference genome fasta file.   
**File output:** Aligned, unsorted BAM file   
 
```
nextflow run bwamem2.nf —params-file <my-params.json> -c <my-nextflow.config> -with-singularity <image.sif>
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
nextflow run sort_and_index.nf —params-file <my-params.json> -c <my-nextflow.config> -with-singularity <image.sif>
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
**File output:** BAM file of marked duplicates and a metrics txt file  
 
```
nextflow run mark_duplicates.nf —params-file <my-params.json> -c <my-nextflow.config> -with-singularity <image.sif>
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

## Variant Calling Workflow and Annotation Workflow [WIP] 

This variant calling worfklow uses Nextflow and is using Singularity. Submodules of the workflow are described below and can be run
independently assuming the necessary input files exist. A parameter file is passed as input using -params-file <my-params.json>, which can begenerated using the templating script. In the example below, parameters are passed as command line argument for to easily demonstrate usage.

### Overview

![variantcall_pipeline](https://github.com/ohsu-cedar-comp-hub/WGS-nextflow-workflow/assets/136844363/1aa03724-3ba7-466d-ab0e-de90ca1d56a3)

### Tools
- Nextflow (DSL2)
- Singularity 
- bcftools/1.10.2
- GATK:4.4.0.0
  - GeneratePileUpSummaries
  - CalculateContamination 
  - Mutect2
  - concat
  - sort
  - index
  - normalize
  - merged
- SNPEff/4.5
- htslib/1.10.2 ## 

### Usage
### Run the variant calling and annotation workflow
```nextflow run WGS-variant-call-annotate-nextflow-pipeline.nf —params-file <my-params.json> -c <my-nextflow.config> -with-singularity <image.sif>```

### 1. GetPileupSummaries
*Generate read counts for estimating of contamination in the tumor tissue using GATK.*  
**File input**: Duplicate-marked BAM files of tumor and matched normal and common germline variants VCF (exac file)  
**File output**: Pileup summary table of ref and alt counts, and allele frequencies at each position

```
nextflow run get_pileup_summaries.nf —params-file <my-params.json> -c <my-nextflow.config> -with-singularity <image.sif>
```

### 2. Calculate Contamination
*Calculate fraction of normal cell contaminants in tumor sample*  
**File input**: Tumor and normal pileup summary tables  
**File output**: Contamination table and tumor segmentation table

```
nextflow run calculate_contamination.nf —params-file <my-params.json> -c <my-nextflow.config> -with-singularity <image.sif>
```


### 3. Run Mutect2 on per chromosome coding sequence files (example given for chrom1 but mutect2 will accept multiple intervals, i.e., all 23 chromosomes as intervals per file)
*Call somatic variants*  
**File input**: Duplicate-marked BAM files of tumor and matched normal, genome reference fasta, germline resource VCF (exac), panel of normals (PoN) VCF, and genomic interval over which to operate  
**File output**: f1r2 .tar.gz file, unfiltered tumor VCF, and VCF stats file


```
nextflow run mutect2.nf —params-file <my-params.json> -c <my-nextflow.config> -with-singularity <image.sif>
```

Script in nextflow:  
```
gatk Mutect2 \
    -R ${mutect_idx} \
    -I ${tumor_bam_sorted} \
    -I ${normal_bam_sorted} \
    --panel-of-normals ${pon} \
    -O ${tumor_bam_sorted.baseName}_unfiltered.vcf \
    --f1r2-tar-gz ${tumor_bam_sorted.baseName}_f1r2.tar.gz \
    -stats ${tumor_bam.baseName}_unfiltered.vcf.stats
```

```docker run quay.io://ohsu-comp-bio/gatk:4.4.0.0 gatk Mutect2 -R reference.fa -I tumor.bam -I normal.bam -normal NORMAL --intervals chr1 -pon gnomad_panel_of_normals.vcf -germline-resource exac_germline_mutation_data.vcf --f1r2-tar-gz "${file%.bam}.f1r2.tar.gz -O "${file%.bam}.unfiltered.vcf"```

### 4. GATK LearnReadOrientationModel
_Learn the read orientation model to refine variant calls by removing technical artifacts._ 
**File input**: f1r2 file from mutect2  
**File output**: Read orientation model .tar.gz file  

```
nextflow run learn_read_orientation.nf —params-file <my-params.json> -c <my-nextflow.config> -with-singularity <image.sif>
```

### 5. Process VCFs 
_Sort, index, normalize and combine (per sample) the VCF files before filtering_`

#### Aggregate across chromosomes with bcftools
**File input**:
**File output**: 

```docker run quay.io//ohsu-comp-bio/bcftools:1.12 bcftools concat -a -f -l listSampleSpecificChromFiles -o "${file%.unf.vcf.gz}.concat.vcf" ```

#####  5a. Sort bgzipped VCFs
``` docker run quay.io//ohsu-comp-bio/bcftools:1.12 bcftools sort "${file}.concat.vcf" -Oz -o  "${file%.unf.concat.vcf.gz}.unfiltered.sorted.vcf.gz"```

###### 5b. Index bgzipped VCF files
``` docker run quay.io//ohsu-comp-bio/bcftools:1.12 bcftools index -t "${file}.unfiltered.sorted.vcf.gz" ```

##### 5c. Normalize pre-filtering for annotation
```docker run quay.io//ohsu-comp-bio/bcftools:1.12 normalize -i "${file}.unfiltered.sorted.vcf.gz" -o "${file%.unfiltered.sorted.vcf.gz}.unfiltered.normalized.vcf"```

#### 6. GATK FilterMutectCalls
_Apply filters to Mutect2 variant calls_ 
```docker run quay.io://ohsu-comp-bio/gatk:4.4.0.0 gatk FilterMutectCalls -R reference.fa -V  "${file}.unfiltered.normalized.vcf" -contamination-table contamination.table --orientation-bias-artifact-priors read-orientation-model.tar.gz -O "${file%.unfiltered.normalized.vcf.gz}.filtered.vcf" ```

#### 7. SNPEff Annotation
_Annotate variants using SNPEff_
```docker run quay.io://ohsu-comp-bio/SNPEff:latest -v genome_version <annotated_variants.vcf> | gzip > annotated_variants.vcf.gz```

