# import pandas as pd
import os
import subprocess

def main():

    # Path to the directory
    fasta_repo = "./fasta_repo"
    out = "./egg_out_3"

    # Loop through files in the directory
    for file in os.listdir(fasta_repo):
        samplename = file[:-6]
        command = f"nice Tools/eggnog-mapper/eggnog-mapper/emapper.py -i {fasta_repo}/{file} -o {out}/{samplename}  --itype CDS  --cpu 14 -m diamond --decorate_gff gtf_file/{samplename}.gtf --decorate_gff_ID_field GeneID --override --tax_scope 33090"
        print(command)
        subprocess.run(command, shell=True, capture_output=True, text=True) #--tax_scope 4072

main()