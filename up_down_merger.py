#!/usr/bin/env python3
"""
Author: Bjorn Wiggers
Description: Up and down differential expression files from 1 input directory. 
Merged into current dir.
Usage: python3 up_down_merger.py *dir*
"""

import pandas as pd
import os
import sys


def main()->None:
    """main function
    Args:
    path: path to the directory you want the script to do things
    Returns:
        None

    """
    samples = []

    path = "." # current dir
    if sys.argv[1]:
        path = sys.argv[1] # unless provided

    for file in os.listdir(path):
        if file.endswith(".tsv"):
            sample_name_a = file.split("_")[0]
            sample_name_b = file.split("_")[1]
            samples.append(sample_name_a + "_" + sample_name_b)
    samples = set(samples)

    for sample in samples:
        print(sample)
        l = []
        a = pd.read_csv(f"{path}/{sample}_downregulatednew.tsv", sep="\t")
        b = pd.read_csv(f"{path}/{sample}_upregulatednew.tsv", sep="\t")
        l.append(a)
        l.append(b)
        result:pd.DataFrame = pd.concat(l)

        result_f = result[["KO_number", "color"]]
        result_f.to_csv(f"{sample}_difeq.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()