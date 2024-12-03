
#import os

# Generate a list of sample names based on the files in the fastqc_test directory
ALL_SAMPLES = [f.replace(".fastq.gz", "") for f in os.listdir("fastqc_test") if f.endswith(".fastq.gz")]
# Extract sample IDs (e.g., SRR24630900, SRR24630901) from the filenames
SAMPLES = glob_wildcards("fastqc_test/{sample}_1.fastq.gz").sample

# This assumes you want to get all .sam files in a folder
SAM_FILES = glob_wildcards("sam_files/{sample}.sam").sample

#Get all .bam files from the folder
BAM_FILES = glob_wildcards("bam_files/{sample}.bam").sample
fastqc_dir = "fastqc_test/fastqc_report"
index_dir = "genome_index"
GENOME = "genome/pepperbase/Capsicum annuum_genome.fasta2.zip"

 
rule all:
	input:
		# expand(fastqc_dir + "/{sample}_fastqc.html", sample=ALL_SAMPLES),
		# "multiqc_report",
		expand("sam_files/{sample}.sam", sample=SAMPLES),
		expand("bam_files/{sample}.bam", sample=SAM_FILES),
		# expand("bam_bai_files/{sample}.bam.bai", sample=BAM_FILES)

# rule fastqc:
# 	input:
# 		fastq = "fastqc_test/{sample}.fastq.gz"
# 	output:
# 		report = fastqc_dir+ "/{sample}_fastqc.html",
# 		zipped_report=fastqc_dir +"/{sample}_fastqc.zip"
# 	shell:
# 		"fastqc {input.fastq} -o {fastqc_dir}"

# rule multiqc:
# 	input:
# 		fastqc_folder="fastqc_test/fastqc_report"
# 	output:
# 		report_dir=directory("multiqc_report")
# 	shell:
# 	        """
# 	        multiqc -f {input.fastqc_folder} -o {output.report_dir}
# 	"""


#rule hisat2_index:
#	input:
#		genome=GENOME
#	output:
#		expand("genome_index/hisat2_Index_base.{i}.ht2", i=range(1, 9))
#	params:
#		threads = 1
#	shell:
#		#"unzip -p '{input.genome}' | hisat2-build -p {params.threads} - genome_index/hisat2_Index_base"
#		"""
#		unzip -p '{input.genome}' > temp_genome.fasta
#		hisat2-build -p 1 temp_genome.fasta genome_index/hisat2_Index_base
#		
#		"""

rule hisat2_mapping:
	input:
		reads_1="data/GSE240943/{sample}_1.fastq.gz",
		reads_2="data/GSE240943/{sample}_2.fastq.gz",
		index="genome/pepperbase/T2T_hisat.1.ht2"
	output:
		"sam_files/{sample}.sam"
	params:
		index="genome/pepperbase/T2T_hisat"
	threads: 1
	shell:
		"hisat2 -p {threads} -x {params.index} -1 {input.reads_1} -2 {input.reads_2} -S {output}"



rule sort_sam: # works for both bowtie2 and hisat2
	input:
		"sam_files/{sample}.sam"
	output:
		"bam_files/{sample}.bam"
	shell:
		"samtools sort {input} -o {output}"



# rule index_bam: # works for both bowtie2 and hisat2
# 	input:
# 		"bam_files/{sample}.bam"
# 	output:
# 		"bam_bai_files/{sample}.bam.bai"
# 	shell:
# 		"samtools index {input} {output}"

#rule convert_gff:
#	input:
#		annotation = ANNOTATION
#	output:
#		temp("genome/pepperbase/Capsicum_annuum_genome.gtf")
#	shell:
#		"gffread {input.annotation} -T -o {output}"

rule stringtie:
	input:
		bamfile="bam_files/{sample}.bam", # stringtie is only used after hisat2
		annotation="genes.gtf"
	output:
		"gtf_file/{sample}.gtf"
	params:
		label="{sample}"                                                                                                                
	shell:
		"stringtie -G {input.annotation} -o {output} -l {params.label} {input.bamfile} -e"


rule prepDE:
    input:
        gtf_file="gtf_file/{sample}"
    output:
        "output/{sample}_gene_count_matrix.csv", 
        "output/{sample}_transcript_count_matrix.csv"
    params:
        script="Tools/prepDE.py",
    shell:
        """
        python {params.script} -i {input.gtf_file}
        """

