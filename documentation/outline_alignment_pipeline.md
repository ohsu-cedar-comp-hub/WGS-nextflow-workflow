## Alignment Pipeline Outline 

Description of the [WGS nextflow alignment pipeline](https://github.com/ohsu-cedar-comp-hub/WGS-nextflow-workflow) with a comparison to the [GDC pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/).

![align_pipeline](https://github.com/ohsu-cedar-comp-hub/WGS-nextflow-workflow/assets/136844363/038d08a6-1621-43e1-b8f9-da62a89871f8)

### Tools 
Versions are just what I have in my conda env as of now. To be updated for consistency with versions in Singularity
- Nextflow (DSL2)
- Singularity
- trimmomatic 0.39
- bwa-mem2 2.2.1
- samtools 1.18
- gatk4 4.3.0.0

|Step | Parameter | COH WGS nextflow pipeline | GDC pipeline |
| --- | --------- | ------------------------- | ------------ |
| 0 | Reference genome | GRCh38.d1.vd1 | GRCh38.d1.vd1 |
| 0 | Read groups | -R flag | Aligned separately then merged << Unsure about this section |
| 1 | BWA algorithm | bwa-mem2 | bwa-mem |
| 2 | BAM sort | samtools | picard SortSam |
| 3 | Mark duplicates | picard MarkDuplicates | picard MarkDuplicates |
| 4 | Co-cleaning workflow | NA | gatk RealignerTargetCreator, IndelRealigner, BaseRecalibrator, PrintReads | 

## Usage  

#### 0. Pre-alignment: quality filtering and trimming  
See https://github.com/ohsu-cedar-comp-hub/WGS-nextflow-workflow/issues/12 << Replace me with real fastqc readme when completed. After QC, trimmomatic is run in paired-end mode.   
**File input**: Paired end reads (read 1 and read 2) and a file of sequencing adapters used in the library prep (TruSeq3-PE.fa)  
**File output**: Paired end reads (trim read 1 and trim read 2), adapters trimmed and quality filtered  
**Call:**  
```
nextflow run trimmomatic.nf —params-file <my-params.json> -c <my-nextflow.config> -with-singularity <image.sif>
```
**Script in nextflow:**
```Shell
trimmomatic \
    PE -phred33 \
    ${read1} \
    ${read2} \
    ${read1.baseName}_1P.fastq.gz \
    ${read1.baseName}_1U.fastq.gz \
    ${read2.baseName}_2P.fastq.gz \
    ${read2.baseName}_2U.fastq.gz \
    ILLUMINACLIP:"${truseq3pefile}":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
```

#### 1. Alignment with Burrows-Wheeler Aligner (BWA) algorithm
 
**File input:** Trimmed, quality filtered, paired-end .fastq files, and an indexed reference genome fasta file.   
**File output:** Aligned, unsorted BAM file   
**Call:**  
```
nextflow run bwamem2.nf —params-file <my-params.json> -c <my-nextflow.config> -with-singularity <image.sif>
```

**Script in nextflow:**
``` Shell
bwa-mem2 mem \
    -K 100000000 -t 6 -Y -M \
    -R "@RG\tID:${params.ID}\tLB:no_library\tPL:illumina\tPU:none\tSM:${trim_read1.baseName}" \
    ${params.idx} ${trim_read1} ${trim_read2} |
    samtools view -Sb -@ 4 > ${trim_read1.baseName}.bam
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

#### 2. Sort and index
**File input**: Unsorted BAM file  
**File output**: Sorted and indexed BAM file (.bam file and .bam.bai file)  
**Call**:

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

#### 3. Mark duplicates
**File input:** Sorted BAM file   
**File output:** BAM file of marked duplicates and a metrics txt file  
**Call**:  
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
