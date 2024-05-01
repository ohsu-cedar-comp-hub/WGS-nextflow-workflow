
## Nextflow WGS tumor evolution Analysis

## Getting Started
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

Test Nexflow successfully installed
```
nextflow info
```

### Example Ceph Bucket Config

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

### Invoking workflow

```
nextflow run hello_world.nf -params-file params.json -c my-config.config
```

Current workflow design for single and paired FASTQ files:
1. QC with Fastqc
2. Trim with TrimmomaticPE
3. Fastqc on trimmed reads
4. MultiQC on all fastqc output files
5. Alignment with bwa-mem2
6. sort and index with samtools
7. Mark duplicates with gatk MarkDuplicates
8. somatic variant calling with gatk Mutect2
9. Learn orientation bias model
10. Get pileup summaries
11. Get segmentation and contamination tables
12. Filter variant calls with gatk Mutect2 (using input files from steps 8-11)
13. Annotate variants with snpEff

