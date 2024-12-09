import pandas as pd
import os
import subprocess

def main():

    # Path to the directory
    fasta_repo = "./fasta_repo"
    out = "./egg_out"

    # Loop through files in the directory
    for file in os.listdir(fasta_repo):
        file_stub = file[:6]
        command = f"nice Tools/eggnog-mapper/emapper.py -i {fasta_repo+"/"+file} -o {out+"/"+file_stub}  --itype CDS  --cpu 8 -m diamond --decorate_gff gtf_file/{file_stub}.gtf --decorate_gff_ID_field GeneID --override"
        subprocess.run(command, shell=True, capture_output=True, text=True)
        # print(file)
    # Print the output

main()