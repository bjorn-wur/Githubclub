
SAMPLES = glob_wildcards("fastqc_test/{sample}.fastq.gz").sample 


# Categorize samples into paired and single-end
paired_samples = list(set(f.split('_')[0] for f in SAMPLES if '_1' in f or '_2' in f))
single_samples = [f for f in SAMPLES if not any(f.startswith(p) for p in paired_samples)]

fastqc_dir = "fastqc_test/fastqc_report"

# Rule to run all jobs
rule all:
	input:
		expand(f"{fastqc_dir}/{{sample}}_fastqc.html", sample=SAMPLES),
		expand(f"{fastqc_dir}/{{sample}}_fastqc.zip", sample=SAMPLES),
		"multiqc_report/multiqc_report.html",
		
		expand("bam_files/{sample}.bam",sample=paired_samples + single_samples),
		expand("gtf_files/{sample}.gtf", sample=paired_samples + single_samples),
		"prepde_input.txt",
		"gene_count_matrix.csv",
		"transcript_count_matrix.csv"

rule fastqc:
        input:
                fastq = "fastqc_test/{sample}.fastq.gz"
        output:
                report = fastqc_dir + "/{sample}_fastqc.html",
                zipped_report = fastqc_dir + "/{sample}_fastqc.zip"
        params:
                fastqc_dir=fastqc_dir
        shell:
                "fastqc {input.fastq} -o {params.fastqc_dir}"

rule multiqc:
	input:
		fastqc_reports = expand("fastqc_test/fastqc_report/{sample}_fastqc.html", sample=SAMPLES),
		fastqc_zips = expand("fastqc_test/fastqc_report/{sample}_fastqc.zip", sample=SAMPLES)
	output:
		"multiqc_report/multiqc_report.html"  # Output path for the multiqc report
	shell:
		"multiqc -f fastqc_test/fastqc_report -o multiqc_report"

rule hisat2_mapping:
	input:
		read_1="fastqc_test/{sample}_1.fastq.gz",
		read_2="fastqc_test/{sample}_2.fastq.gz"
	output:
		temp("{sample}.sam")
	params:
		index="genome/pepperbase/T2T_hisat"
	threads: 4
	shell:
		"""
		hisat2 -p {threads} -x {params.index} -1 {input.read_1} -2 {input.read_2} -S {output}
		"""

# Rule for single-end mapping
rule hisat2_single_end_mapping:
	input:
		"fastqc_test/{sample}.fastq.gz"
	output:
		temp("{sample}.sam")
	params:
		index="genome/pepperbase/T2T_hisat"
	threads: 4
	shell:
		"""
		hisat2 -p {threads} -x {params.index} -U {input} -S {output}
		"""

rule sort_sam: # works for both bowtie2 and hisat2
	input:
		"{sample}.sam"
	output:
		"bam_files/{sample}.bam"
	shell:
		"samtools sort -o {output} {input}"


rule index_bam: # works for both bowtie2 and hisat2
	input:
		"bam_files/{sample}.bam"
	output:
		"bam_files/{sample}.bam.bai"
	shell:
		"samtools index {input}"

rule stringtie:
	input:
		bamfile="bam_files/{sample}.bam", # stringtie is only used after hisat2
		annotation="Capsicum_annuum_genome_fixed.gtf"
	output:
		"gtf_files/{sample}.gtf"
	params:
		label="{sample}"                                                                                                                
	shell:
		"stringtie -G {input.annotation} -o {output} -l {params.label} {input.bamfile} -e"


rule prepde_input:
	output:
		a="prepde_input.txt"
	run:
        # Retrieve all sample names from BAM files
		BAM_FILES = glob_wildcards("bam_files/{sample}.bam").sample
        # Open the output file
		with open(output[0], "w") as f:
			for i, sample in enumerate(BAM_FILES):
				f.write(f"{sample} gtf_files/{sample}.gtf\n")

rule prepDE:
	input:
		prepde="prepde_input.txt",
	output:
		"gene_count_matrix.csv", 
		"transcript_count_matrix.csv"
	params:
		script="Tools/prepDE.py"
	shell:
		"""
		python {params.script} -i {input.prepde}
		"""

