#!/usr/bin/env python3
"""
Author: Ludovica Piccioli
Description: merges annotation files into one file, parse TSV files and creates
    a new TSV file each with the relevant information needed
Usage: python3 parse_tsv.py
"""

from pathlib import Path

def merge_annotation_files(list_files):
    """Merges annotation files into one file

        Args:
            list_files (list): a list of annotation files

        Returns:
            None

    """
    list_geneids = []
    for i in range(len(list_files)):
        with open(list_files[i], "r") as fin:
            # skip the first and last 4 lines, keep headers
            lines = fin.readlines()[4:-4]
            if i == 0: # if it is the fist annotation file
                for j in range(len(lines)):
                    gene_id = lines[j].split("\t")[11] # gene id is 12th column
                    list_geneids.append(gene_id)
                with open("annotation_merged", "w") as fout:
                    fout.writelines(lines) #write the annotations
            else: # if it is not the first file
                with open("annotation_merged", "a") as fout:
                    for j in range(len(lines)):
                        gene_id = lines[j].split("\t")[11]
                        #if the annotation for that gene is not already present
                        if gene_id not in list_geneids:
                            list_geneids.append(gene_id)
                            fout.write(lines[j])

def get_color(log_fold_change):
    """Get a color based on the log fold 2 change

        Args:
            log_fold_change (float): the log fold 2 change

        Returns:
            color (str): the assigned color

    """
    list_colors_up = ["#e8aa9f", "#e29182", "#db7765", "#d45d48", "#ca462e",
                      "#ad3c28","#903221", "#73281a", "#561e14", "#39140d"]
    list_colors_down = ["#96ddf2", "#76d3ee", "#56c9ea", "#36bee6", "#1ab2df",
                        "#1698bf", "#137f9f", "#0f657f", "#0b4c5f", "#07323f"]
    color = None
    for n in range(10):
        if 0+n <= log_fold_change  < 1+n:
            color = list_colors_up[0+n]
        elif 0-n > log_fold_change >= -1-n:
            color = list_colors_down[0+n]
    return color


def get_ko_numbers(annotation_file):
    """Get the KO numbers for a gene id from the annotation file

            Args:
            annotation file (str): the name of the annotation file

            Returns:
                id_ko (dict): dictionary with gene ids and corresponding ko

        """
    with open(annotation_file, "r") as fin:
        lines = fin.readlines()
    #ko numbers are in column index 11
    id_ko = {}
    for item in lines:
        line = item.split("\t")
        # do not include '-mrna' in the gene id
        gene_id, ko = line[0][:-5], line[11]
        # do not include genes for which the ko is not provided
        if gene_id not in id_ko and ko != "-":
            # if there are multiple ko numbers for the gene
            if "," in ko:
                ko_multiple = ko.split(",")
                kos = []
                # only extract the number itself
                for number in ko_multiple:
                    ko_num = number[3:]
                    kos.append(ko_num)
                id_ko[gene_id] = kos
            elif "," not in ko:
                id_ko[gene_id] = ko[3:]
    return id_ko


def parse_tsv(tsv_file, dict_ids_kos, new_tsv):
    """Parses a tsv file and writes a new tsv

            Args:
            tsv (str): the name of the tsv file
            dict_ids_kos (dict): a dictionary with gene ids and ko numbers
            new_tsv (str): the name of the new tsv file

            Returns:
                None

        """
    with open(tsv_file, "r") as fin:
        # skip the header line
        lines = fin.readlines()[1:]
    with open(new_tsv, "w") as fout:
        fout.write("GeneID\tlog2FoldChange\tcolor\tKO_number\n")
        for item in lines:
            line = item.split("\t")
            gene_id = line[0].strip('"')
            log_fold_change = line[2]
            if gene_id in dict_ids_kos:
                ko_n = dict_ids_kos[gene_id]
                color = get_color(float(log_fold_change))
                fout.write(f"{gene_id}\t{float(log_fold_change):.3f}\t"
                           f"{color}\t")
                if isinstance(ko_n, list): # if there are multiple kos
                    quan = len(ko_n)
                    for i in range(quan):
                        if i != quan-1:
                            fout.write(ko_n[i]+ ", ")
                        else:
                            fout.write(ko_n[i])
                elif isinstance(ko_n, str):
                    fout.write(ko_n)
                fout.write("\n")


def main():
    """This is the main function
    """
    annotation_files = ["SRR24630893.emapper.annotations",
                        "SRR24630912.emapper.annotations",
                        "SRR24630931.emapper.annotations",
                        "SRR24630932.emapper.annotations",
                        "SRR24630948.emapper.annotations",
                        "SRR24630949.emapper.annotations",
                        "SRR24630967.emapper.annotations",
                        "SRR25031739.emapper.annotations"]
    list_files = ["p_capsici_4h_downregulated.tsv","TMV_4h_downregulated.tsv",
                  "Xag8ra_3h_downregulated.tsv", "XCV1_3h_downregulated.tsv",
                  "XCV3_3h_downregulated.tsv", "p_capsici_4h_upregulated.tsv",
                  "TMV_4h_upregulated.tsv", "Xag8ra_3h_upregulated.tsv",
                  "XCV1_3h_upregulated.tsv", "XCV3_3h_upregulated.tsv"]
    merge_annotation_files(annotation_files)
    dict_kos = get_ko_numbers("annotation_merged")
    for file in list_files:
        filename = "R_output/tsv_output/" + file
        path_obj = Path(filename)
        new_tsv = "tryingoutput/" + path_obj.stem + "new.tsv"
        parse_tsv(filename, dict_kos, new_tsv)

if __name__ == "__main__":
    main()