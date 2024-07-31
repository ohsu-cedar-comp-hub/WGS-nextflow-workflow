# WGS Nextflow Workflow

This Nextflow workflow aligns, performs somatic variant calls, and annotates variants on whole genome sequencing matched tumor-normal data. The pipeline is intended to run with a container using Singularity. Submodules of the workflow are described below.

## Getting Started

Check that Java version is 11 through 21. Nextflow has been checked to run with Java 11-21. Otherwise install Java 17 via SDKMAN shown [here](https://www.nextflow.io/docs/latest/install.html).  

```
java -version
```

### Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) 

Option 1: Install with conda/mamba/micromamba (recommended)

```
# create a conda environment
conda create -n my_nextflow_environment -c bioconda nextflow

# activate conda environemnt
conda activate my_nextflow_environment

# confirm install
nextflow info

# update
nextflow self-update
```

Option 2: Install with curl (not for use on Exacloud)

```
# Move to a directory you want Nextflow to sit on your system
cd <your-directory-of-choice>

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

### Pull Singularity images
[Read more about Singularity here.](https://docs.sylabs.io/guides/3.2/user-guide/cli/singularity_pull.html)
 
Dockerfiles for building images are provided for each tool in the pipeline in the tool directory, `tools/<tool name>/Dockerfile`. Docker images have been pre-built and are hosted on quay.io at the following links:

- fastqc https://quay.io/repository/ohsu-comp-bio/fastqc
- multiqc https://quay.io/repository/ohsu-comp-bio/multiqc
- trimmomatic https://quay.io/repository/ohsu-comp-bio/trimmomatic
- bwa https://quay.io/repository/ohsu-comp-bio/bwa 
- samtools https://quay.io/repository/ohsu-comp-bio/samtools
- gatk https://quay.io/repository/ohsu-comp-bio/gatk
- bcftools https://quay.io/repository/ohsu-comp-bio/bcftools
- snpeff https://quay.io/repository/ohsu-comp-bio/snpeff

Because Docker cannot be easily used on Exacloud, images can be pulled using Singularity, which converts them to .sif images. To pull images:

```
## Load singularity module 
module load /etc/modulefiles/singularity/current

## set singularity cache directory variable; can live anywhere but gscratch is recommended.
export SINGULARITY_CACHEDIR=/home/exacloud/gscratch/CEDAR/<user>/singularity/cache

## set temporary directory variable
export TMPDIR=<your current working directory>/tmp

## Navigate to the directory you would like your .sif images to be
cd <sif file dir>

## Pull your images
singularity pull <name_of_singularity_image>.sif docker://quay.io/ohsu-comp-bio/<name of tool>
```

Alternatively, images for this pipeline are temporarily hosted at `/home/groups/CEDAR/lancasru/WGS_COH_NF/config_sif`.

#### Set up parameters file

See an example template here: [params](https://github.com/ohsu-cedar-comp-hub/WGS-nextflow-workflow/blob/2d5bca446ab3619a481a9f3c24115f708cd249cc/example/nf_channels_params_template.json).

The contents of the params file are shown below; add your customized paths to the following variables:
```Json
{
"container_fastqc" : "/path/to/singularity/image/fastqc.sif",
"container_multiqc" : "/path/to/singularity/image/multiqc.sif",
"container_trimmomatic" : "/path/to/singularity/image/trimmomatic.sif",
"container_samtools" : "/path/to/singularity/image/samtools.sif",
"container_gatk" : "/path/to/singularity/image/gatk.sif",
"container_bwa" : "/path/to/singularity/image/bwa.sif",
"container_snpeff" : "/path/to/singularity/image/snpeff.sif",
"container_bcftools" : "/path/to/singularity/image/bcftools.sif",

"all_read_pairs" : "/path/to/raw/fastq/*_R{1,2}*.fastq.gz",
"id" : "UUID-UUID-UUID",
"outdir" : "/path/to/output/directory",

"truseq3pefile" : "path/to/adapters/fasta",
"idx" : "/path/to/reference/index/for/alignment",
"exac" : "/path/to/common/germline/variants/vcf",
"gnomad" : "/path/to/germline/resource/gnomad",
"pon" : "/path/to/panel/of/normals",
"mutect_idx" : "/path/to/indexed/and/dictionaried/reference",
}
```

**all_read_pairs**  
The `all_read_pairs` variable points to your raw fastq files for a single patient. The pipeline expects an input of 4 files with the sample ID separated from the rest of the file name by an underscore `_`, as below:

```
├── my_fastq_files
    ├── YOURSAMPLEID_G_R1.fastq.gz
    ├── YOURSAMPLEID_G_R2.fastq.gz
    ├── YOURSAMPLEID_T_R1.fastq.gz
    ├── YOURSAMPLEID_T_R2.fastq.gz
```
_Filename requirements_

Regular expressions expect filenames for raw fastq files to be structured like this:
```
SAMPLEID_{G,T}_R{1,2}_additionalmetainfo.fastq.gz
```

Where {G,T} is either G (germline) or T (tumor), and {1,2} is either read 1 or 2. This naming scheme works with the suggested regular expression in the `all_read_pairs` variable; for other file names, edit the `*_R{1,2}.fastq.gz` regular expression in the params file to encompass your naming scheme. This stringent naming requirement will be updated in future releases. Currently, multiple tumors per patient are not supported. 

**id**  
A unique identifier to be included in your read group header

**outdir**  
This is your output directory. Make this empty directory before running the pipeline.

**truseq3pefile**  
A FASTA file of adapters used in sequencing. This is the file referenced by trimmomatic when trimming adapters.

**idx**  
The path to your reference genome index; eg, when using GRCh38.d1.vd1.fa, the path will be `path/to/reference/GRCh38.d1.vd1`. The index is generated by running `bwa index <reference genome>.fa`. The output includes index files with suffixes like .amb, .ann, .bwt, .pac, etc. Read more about indexing the reference genome [here](https://bio-bwa.sourceforge.net/bwa.shtml). 

**exac**  
Path to common germline variants VCF file. GATK best practices recommends using either gs://gatk-best-practices/somatic-b37/small_exac_common_3.vcf or gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz, depending on your reference, from their [resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle), available [on Google Cloud](https://console.cloud.google.com/storage/browser/gatk-best-practices/). 

**gnomad**  
Path to germline resource; when working with human data this is almost guarunteed to be from gnomad. GATK best practices recommends using af-only-gnomad.hg38.vcf.gz, which is [a copy of the gnomAD VCF stripped of all unnecessary INFO fields](https://gatk.broadinstitute.org/hc/en-us/articles/360050722212-FAQ-for-Mutect2). It is also available [on Google Cloud](https://console.cloud.google.com/storage/browser/gatk-best-practices/). 

**pon**  
Path to panel of normals. GATK best pratices recommends using gatk4_mutect2_4136_pon.vcf, available [on Google Cloud](https://console.cloud.google.com/storage/browser/gatk-best-practices/).  

**mutect_idx**  
Path to reference genome for running Mutect2 (includes the reference.fa, reference.dict, and reference.fa.fai). This is different from your reference genome specified for bwa-mem2 alignment, because the indices differ. This is available from [Google Cloud](https://console.cloud.google.com/storage/browser/gatk-best-practices/) or by running `samtools faidx reference.fa`


### Set up config file

User specifications can be referenced in config files, and Nextflow will pull information from them in order of preference described [here](https://www.nextflow.io/docs/latest/config.html). For specificity, we suggest running Nextflow with the command `-c path/to/nextflow.config` to ensure you are referencing the correct config file. 

Below is an example of a config file that is set up to run in Singularity containers and use the slurm executor. Memory and time should be adjusted as necessary. 

```Nextflow
singularity {
    enabled = true
    autoMounts = true
    cacheDir = '/home/exacloud/gscratch/CEDAR/lancasru/singularity'
    }

// specify executor args for each process

process {
    executor = 'slurm'
    queue = 'exacloud'

    withName: 'MARKDUPLICATES' {
        memory = '40 GB'
        time = '01:00:00'
    }
    withName: 'TRIMMOMATICPE' {
        memory = '20 GB'
        time = '00:30:00'
    }
    withName: 'FASTQCRAW' {
        memory = '2 GB'
        time = '00:10:00'
    }
    withName: 'FASTQCTRIM' {
        memory = '2 GB'
        time = '00:10:00'
    }
    withName: 'BWAMEM2' {
        memory = '42 GB'
        time = '01:00:00'
        cpus = 2
    }
    withName: 'MULTIQC' {
        memory = '2 GB'
        time = '00:10:00'
    }
    withName: 'SORT' {
        memory = '8 GB'
        time = '00:15:00'
    }
    withName: 'SORTANDINDEX' {
        memory = '8 GB'
        time = '00:15:00'
    }
    withName: 'GETPILEUPSUMMARIES' {
        memory = '2 GB'
        time = '00:15:00'
    }
    withName: 'CALCULATECONTAMINATION' {
        memory = '2 GB'
        time = '00:15:00'
    }
    withName: 'MUTECT2' {
        memory = '40 GB'
        time = '00:45:00'
    }
    withName: 'BGZIP' { 
        memory = '2 GB'
        time = '00:15:00'
    }
    withName: 'PREPAREVCF' {
        memory = '2 GB'
        time = '00:15:00'
    }
    withName: 'MERGESTATS' {
        memory = '2 GB'
        time = '00:15:00'
    }
    withName: 'LEARNORIENTATION' {
        memory = '2 GB'
        time = '00:15:00'
    }
    withName: 'FILTERMUTECT' {
        memory = '40 GB'
        time = '00:45:00'
    }
    withName: 'ANNOTATE' {
        memory = '20 GB'
        time = '00:45:00'
    }
}
```

Ensure that you have specified all the individual memory/time requirements in your config file for each process. These will vary based on size of the fastq file. 

### Set up slurm job script

Nextflow can use the slurm executor to submit jobs; one job per process. To run jobs on Exacloud, submit a job script (sbatch). 

When running Nextflow, you will need:
- the nextflow script `(workflows/<step>/<script name>.nf)`
- the parameters file
- the nextflow config file 

Make sure all of these args are defined in your job script.  
As long as memory has been specified for each process, the sbatch directive in the job script only requires the memory needed for Nextflow run to submit your jobs to slurm. In this example, it is 1G; however, it could be lower. Some documentation suggests allocating at least 128MB. Make sure your allocated time is sufficient as this job needs to be running the whole time in order for Nextflow to submit all of the processes in the workflow. Depending on the size of your files, this will easily be upwards of 12 hours. For jobs that exceed 36 hours, add the `##SBATCH qos=long_jobs` directive.

```Shell
#!/bin/bash

#SBATCH --mem=1G
#SBATCH --time=08:30:00    ## time=<dd:hh:mm:ss> 
#SBATCH --job-name="my_nextflow_job"
#SBATCH --partition "exacloud"

nextflow run path/to/WGS-nextflow-workflow/workflows/alignment/align.nf -c path/to/nextflow.config -params-file path/to/params-file.json
```

There are 15 jobs spawned by the alignment workflow for one sample:
```
FASTQCRAW (2)
TRIMMOMATICPE (2)
FASTQCTRIM (2)
MULTIQC
BWAMEM2 (2)
SORT (2)
MARKDUPLICATES (2)
SORTANDINDEX (2)
```

There are 53 jobs spawned by the somatic variant calling workflow for one sample: 
```
MUTECT2 (23)
BGZIP (23)
GETPILEUPSUMMARIES
CALCULATECONTAMINATION
MERGESTATS
LEARNORIENTATION
PREPAREVCF
FILTERMUTECT
ANNOTATE
```

Given the number of jobs spawned, it may be in your interest to limit the amount of jobs submitted in parallel by Nextflow by setting `executor.queueSize` in the config file. There are several more options you can set detailed [here](https://www.nextflow.io/docs/latest/config.html#scope-executor).

If your job is killed, it is easy to resume your job from where you left off by adding the `resume <session ID>` argument to the command, eg., 

```
nextflow run path/to/WGS-nextflow-workflow/workflows/alignment/align.nf -c path/to/nextflow.config -params-file path/to/params-file.json -resume 4dc656d2-c410-44c8-bc32-7dd0ea87bebf
```

The session ID is a unique string of hyphenated numbers and characters that can be found in your nextflow.log file. 


## Running the workflow

Activate your Nextflow conda environment, load the singularity module, and submit your slurm job script. It is expected that warnings will be produced for processes that are specified in your config file but not actively used in the workflow (ie, when running alignment, processes for variant calling will generate warnings and vice versa). These warnings can safely be ignored. 

```
# Activate nextflow conda env
conda activate my_nextflow_environment

# move to directory where you want your nextflow work directory and logs to be created
cd <my nextflow run>

# Load singularity module
module load /etc/modulefiles/singularity/current

# submit your slurm job script
sbatch my_nextflow_script.srun
```

---

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

![variantcall_pipeline](https://source.ohsu.edu/storage/user/478/files/60f79880-9bb4-4543-960c-87bcd13ce1a6)

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

### 7. snpEff Annotation
_Annotate variants using SNPEff_  
**File input**: Filtered VCF and genome version GRCh38.86 from snpEff pre-built database  
**File output**: Annotated variants VCF

Script in nextflow:

```
java -Xmx8g -jar /usr/src/app/snpEff/snpEff.jar GRCh38.86 ${filtered_vcf} -cancer > ${sample_id}_annotated.vcf
```

### 8. SnpSift filtering
_Filter out variants not marked as "PASS" by FilterMutectCalls and under the specified allelic depth_
**File input**: Annotated variants VCF
**File output**: Filtered, annotated VCF file

Script in nextflow:

```
cat ${filtered_vcf} | java -Xmx8g -jar /usr/src/app/snpEff/SnpSift.jar filter "( ((GEN[0].AD[0] >= 3) | (GEN[0].AD[1] >= 4)) | ((GEN[1].AD[0] >= 3) | (GEN[1].AD[1] >= 4)) ) & ( ( na FILTER ) | (FILTER = 'PASS') )" > ${sample_id}_annotated_PASSED.vcf
```

### Output directory structure 

Running the pipeline generates a total of 19 files, structured like this:
```
aligned/
	|- unsorted
			|- SAMPLE_G.bam
			|- SAMPLE_T.bam
	|- duplicate_marked
			|- SAMPLE_G_sorted_marked_duplicates_sorted.bam
			|- SAMPLE_G_sorted_marked_duplicates_sorted.bam.bai
			|- SAMPLE_T_sorted_marked_duplicates_sorted.bam
			|- SAMPLE_T_sorted_marked_duplicates_sorted.bam.bai
annotated/
	|- SAMPLE_T_annotated.vcf
	|- SAMPLE_T_annotated_PASSED.vcf
multiqc/
	|- multiqc_report.html
summaries/
	|- SAMPLE_G_sorted_marked_duplicates_sorted.getpileupsummaries.table
	|- SAMPLE_T_sorted_marked_duplicates_sorted.getpileupsummaries.table
svc/
	|- SAMPLE_T_filtered.vcf
	|- SAMPLE_T_read-orientation-model.tar.gz
	|- SAMPLE_T_unfiltered.vcf.all.stats
	|- SAMPLE_T_unfiltered.vcf.gz
	|- SAMPLE_T_unfiltered_normalized_sorted.vcf.gz
	|- SAMPLE_T_unfiltered_normalized_sorted.vcf.gz.tbi
tables/
	|- SAMPLE_T_sorted_marked_duplicates_sorted_contamination_table
	|- SAMPLE_T_sorted_marked_duplicates_sorted_segmentation_table
```
