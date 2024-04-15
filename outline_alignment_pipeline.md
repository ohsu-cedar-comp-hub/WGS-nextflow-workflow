## Alignment Pipeline Outline 

Description of the [COH WGS nextflow pipeline](https://github.com/ohsu-cedar-comp-hub/WGS-nextflow-workflow) in comparison to the [GDC pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/).

| Parameter | COH WGS nextflow pipeline | GDC pipeline |
| --------- | ------------------------- | ------------ |
| File input | FASTQ | FASTQ or BAM |
| File output | BAM | BAM |
| Reference genome | GRCh38.d1.vd1 | GRCh38.d1.vd1 |
| Read groups | -R flag | Aligned separately then merged << Unsure about this section |
| BWA algorithm (1) | bwa-mem2 | bwa-mem |
| BAM sort (2) | samtools | picard SortSam |
| Mark duplicates | picard MarkDuplicates | picard MarkDuplicates | 

#### 0. Pre-Alignment: see FastQC workflow
https://github.com/ohsu-cedar-comp-hub/WGS-nextflow-workflow/issues/12 << Replace me with real readme when completed

#### 1. Alignment with Burrows-Wheeler Aligner (BWA) algorithms
BWA-MEM is used if mean read length is greater than or equal to 70 bp.    
**File input:** Trimmed, quality filtered, paired-end .fastq files (assigned to variable "read 1" and "read 2" in parameters file).    
**File output:** Aligned .bam file   
**Call:** 
``` Shell
bwa-mem2 mem \
-K 100000000 -t 6 -Y -M \
-R "@RG\tID:${params.ID}\tLB:no_library\tPL:illumina\tPU:none\tSM:${trim_read1.baseName}" \
${params.idx} ${trim_read1} ${trim_read2} |
samtools view -Sb -@ 4 > ${trim_read1.baseName}.bam
```
**Comparison GDC call**:  
```Shell
bwa mem \
-t 8 \
-T 0 \
-R <read_group> \
<reference> \
<fastq_1.fq.gz> \
<fastq_2.fq.gz> |
samtools view \
-Shb
-o <output.bam> -
```

#### 2. Sort and index
**File input**: Unsorted BAM file  
**File output**: Sorted and indexed BAM file (.bam file and .bam.bai file)

#### 3. Mark duplicates
