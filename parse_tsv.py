#!/usr/bin/env python3
"""
Author: Ludovica Piccioli
Description: parses a FASTA file to only keep the reads present in the TSV file
Usage: python3 parse_tsv.py <tsv_up> <tsv_down> <fasta_file> <new_fasta>
"""
import os
from sys import argv

def parse_tsv(tsv_file, gene_ids):
    """Parses a TSV file and makes a list of gene ids

        Args:
            tsv_file (string): input TSV file
            gene_ids (list): a list of gene IDs

        Returns:
            gene_ids (list): the expanded list of gene IDs

    """
    with open(tsv_file, "r") as fin:
        lines = fin.readlines()[1:] #skip the header line
    for i in range(len(lines)):
        line = lines[i].split("\t")
        gene_id = line[0].strip('"')
        if gene_id not in gene_ids:
            gene_ids.append(gene_id)
    return gene_ids


def parse_fasta(fasta_file, gene_ids, new_fasta):
    """Parses a FASTA file and writes a new file with the sequences of interest

            Args:
                fasta_file (string): input FASTSA file
                gene_ids (list): a list of gene IDs of interest
                new_fasta (string): the name of the new fasta file

            Returns:
                None

        """
    with open(fasta_file, "r") as fin:
        lines = fin.readlines()

    with open(new_fasta, "w") as fout:
        for i in range(len(lines)):
            if lines[i].startswith(">"):
                # print(lines[i])
                line = lines[i].strip().split("-")[0]
                gene_id = line[1:]
                # print(gene_id)
                if gene_id in gene_ids:
                    fout.write(lines[i])
                    for n in range(1, len(lines) - i):
                        if not lines[i + n].startswith(">"):
                            fout.write(lines[i + n])
                        else:
                            break


def main():
    """This is the main function
    """
        # Path to the directory
    r_output_dir = "./fasta_repo"
    tsv_up, tsv_down, fasta_file, new_fasta = argv[1], argv[2], argv[3], argv[4]
    tsv_up = "up_reg.tsv"
    tsv_down = "down_reg.tsv"
    gene_ids_empty = []
    # Loop through files in the directory
    for file in os.listdir(r_output_dir):
        file2 = "./fasta_repo/" + file
        fasta_file = file2
        new_fasta = f"{file}"
        # command = f"python3 parse_tsv.py up_reg.tsv down_reg.tsv {file} test_48/{file}"
        # subprocess.run(command, shell=True, capture_output=True, text=True)
        gene_ids_up = parse_tsv(tsv_up, gene_ids_empty)
        gene_ids_tot = parse_tsv(tsv_down, gene_ids_up)
        parse_fasta(fasta_file, gene_ids_tot, new_fasta)




    # gene_ids_tot = parse_tsv(tsv_down, gene_ids_up)
    # parse_fasta(fasta_file, gene_ids_tot, new_fasta)


if __name__ == "__main__":
    main()
