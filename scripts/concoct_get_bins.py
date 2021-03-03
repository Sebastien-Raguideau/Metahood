#!/usr/bin/env python3

from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
import argparse
import numpy as np
from collections import Counter

def get_initial_number_of_bins(file, MAX_BIN_NB):
    nb_bin=int(2*np.median(list(Counter([header.split(" ")[1] for header,seq in sfp(open(file))]).values())))
    return min(nb_bin,int(MAX_BIN_NB))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", help="SCG file")
    parser.add_argument("-m", help="max_bin_nb")
    args = parser.parse_args()
    Fna_file=args.s
    Max_bin=args.m

    min_val=get_initial_number_of_bins(Fna_file, Max_bin)
    print(min_val)
