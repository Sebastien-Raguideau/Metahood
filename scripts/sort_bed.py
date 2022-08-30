#!/usr/bin/env python
# -*- coding: latin-1 -*-

from collections import defaultdict,Counter
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
import numpy as np
import argparse


def main(contig_file,bed_file,out_file,genome_file):
    index_contig = {}
    contig_index = {}
    for index,(header,seq) in enumerate(sfp(open(contig_file))):
        contig = header.split(" ")[0]
        contig_index[contig] = [index,len(seq)]
        index_contig[index] = contig
    bed_lines = [line.rstrip() for line in open(bed_file)]
    new_lines =  sorted(bed_lines,key=lambda x:contig_index[x.split("\t")[0]][0])
    with open(out_file,"w") as handle:
        handle.writelines("%s\n"%line for line in new_lines)
    if genome_file:
        with open(genome_file,"w") as handle:
            handle.writelines("%s\t%s\n"%(index_contig[index],contig_index[index_contig[index]][1]) for index in range(len(index_contig)))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bed_file", help="bed file to sort")
    parser.add_argument("contigs", help="contig file used for order")
    parser.add_argument("out", help="new bed file name")
    parser.add_argument("-g", help="create a genome file for bamtools",default = "")
    args = parser.parse_args()
    main(args.contigs,args.bed_file,args.out,args.g)
