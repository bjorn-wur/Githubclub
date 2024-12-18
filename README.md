# BIF Template

Welcome to BIF template, please follow the steps below.

## Project Description

Template for developing Python applications.

## Project Information

- **Name**: template thingy 
- **Description**: Simple python project template.
- **Version**: 0.0.1
- **Author**: Bjorn Wiggers

## Quickstart guide:

## How to run:
The snakemake pipeline can be used to run the analysis up until the DESEQ2 which is done through an R script, which can be run as follows:

(Assuming you have a gff file and a fasta file containg the whole genome.)
From this,.we require a complete list of all th genes present in the original genome, this way it can be annotated using EGGNOG. If you want to regenerate this file you would use the following command:

Tools/gffread/gffread Capsicum_annuum_genome.gff -g "Capsicum_annuum_genome.fasta2" -w total_genome.fasta

EGGNOG IS NOT PRE-INSTALLED IN THIS REPO, IT HAS TO BE ADDED MANUALLY DUE TO SIZE.
-- eggnog can be downloaded from the following link, and should be installedin the ./Tools directory 
    (so the full path to the final script: Tools/eggnog-mapper/eggnog-mapper/emapper.py)
https://github.com/eggnogdb/eggnog-mapper

-- the following script from the eggnog repo has to be run before running eggnog itself
download_eggnog_data.py

In order to get associated data with each gen eggnog is used on the resulting total_genome fasta file with the following command: 
Tools/eggnog-mapper/eggnogmapper/emapper.py -i total_genome.fasta -o eggnog_output  --itype CDS  --cpu 4 -m diamond --decorate_gff Capsicum_annuum_genome.gtf --decorate_gff_ID_field GeneID --override --tax_scope 33090

Using the annotation_and_tsv_parsing.py this would then create a new file containing the:
GeneID	log2FoldChange	color	KO_number 
p_capsici_4h_downregulatednew.tsv
(usage python3 annotation_and_tsv_parsing.py folder_that_contains_r_output)

this finally can be combined for the pathway review on kegg using the up_down_merger.py, to create a singular file for both up and down regulated samples. 


