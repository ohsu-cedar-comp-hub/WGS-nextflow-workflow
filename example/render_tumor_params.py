import jinja2

template_path = './tumor_params_file_template.json'

params = {
    "read1": "path/to/fastq_R1.fastq.gz",
    "read2": "path/to/fastq_R2.fastq.gz",
    "fastqc_read1": "path/to/fastqc_R1_report/fastqc.zip",
    "fastqc_read2": "path/to/fastqc_R2_report/fastqc.zip",
    "trim_read1": "path/to/trim_read1/1P.fastq.gz",
    "trim_read2": "path/to/trim_read2/2P.fastq.gz",
    "trim_fastqc_read1": "path/to/fastqc_trim_read1_report/1P_fastqc.zip",
    "trim_fastqc_read2": "path/to/fastqc_trim_read2_report/2P_fastqc.zip",
    "outdir": "path/to/output/directory",
    "truseq3pefile": "path/to/TruSeq3-PE.fa",
    "idx": "path/to/bwa_mem2_idx/GRCh38.bm2.d1",
    "pon": "path/to/pon.vcf",
    "id": "unique_str_id",
    "normal_bam_unsorted": "path/to/normal_unsorted_bam.bam",
    "normal_bam_sorted": "path/to/normal_sorted_bam.bam",
    "exac": "path/to/exac.exac"
}

env = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath='.'))
template = env.get_template(template_path)

output = template.render(params)

with open('tumor_params_file.json', 'w') as f:
    f.write(output)

print("Tumor params file created successfully.")

