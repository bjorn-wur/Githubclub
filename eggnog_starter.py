#!/usr/bin/env python3
"""
Author: Bjorn Wiggers
Description: loops over input fasta file and eggnogs maps it to a ./egg_out directory directory.
Usage: python3 eggnog_starter.py *dir*

please note the following requirements:
    output directory has to exist
    input directory has to exist. 
    a directory of gtf files in 
    gtf_file/{samplename}.gtf 
    has to exist with the same basic file names as the fasta files.
    nice command has to exist..
"""

import os
import subprocess
import sys

def main():

    # Path to the directory
    fasta_repo = "./fasta_repo"
    out = "./egg_out"

    # Loop through files in the directory
    for file in os.listdir(fasta_repo):
        samplename = file[:-6]
        command = f"nice Tools/eggnog-mapper/eggnog-mapper/emapper.py -i {fasta_repo}/{file} -o {out}/{samplename}  --itype CDS  --cpu 4 -m diamond --decorate_gff gtf_file/{samplename}.gtf --decorate_gff_ID_field GeneID --override --tax_scope 33090"
        print(command)
        subprocess.run(command, shell=True, capture_output=True, text=True) #--tax_scope 4072

main()