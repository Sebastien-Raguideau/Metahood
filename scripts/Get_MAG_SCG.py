#!/usr/bin/env python

import argparse
import sys
import glob
import os

from collections import defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("cog_annotation", help="annotation of the full assembly")
    parser.add_argument("cluster_definition", help="csv file, firts column is the contig, second is the bin")
    parser.add_argument("mag_list", help="list of mags")    
    parser.add_argument("cogs_list", help="list of COGs to collate from bins")
    parser.add_argument(
        "output_folder", help="output folder name, where do we want to put the sequences")
    args = parser.parse_args()
    # set of SCG
    with open(args.cogs_list, 'r') as f:
        core_cogs = set([x.rstrip() for x in f.readlines()])
    # get mags 
    mags={element.rstrip() for element in open(args.mag_list)}
    # get mags definition
    contigs_to_mags={}
    with open(args.cluster_definition) as handle:
        for line in handle:
            contig,bin_=line.rstrip().split(",")
            if bin_ in mags:
                contigs_to_mags[contig]=bin_
    # get SCGs linked to bin
    corecog_bin = defaultdict(lambda:defaultdict(lambda:["",""]))
    for header,seq in sfp(open(args.cog_annotation)):
        orf,COG,strand=header.rstrip().split(' ')
        contig="_".join(orf.split("_")[:-1])
        if contig in contigs_to_mags:
            if len(seq)>len(corecog_bin[COG][contigs_to_mags[contig]][1]) :
                corecog_bin[COG][contigs_to_mags[contig]]=[header,seq]
    # output
    path_output = args.output_folder
    if not os.path.isdir(path_output):
        os.system("mkdir -p "+path_output)
    for corecog, dict_List_bins in corecog_bin.items():
        with open(path_output + corecog + ".fna", "w") as output_handle:
            for bin_name, (header, seq) in dict_List_bins.items():
                new_name = bin_name + "_" + corecog+" "+header
                output_handle.write(">"+new_name+"\n"+seq+"\n")


if __name__ == "__main__":
    # TODO consider parsing arguments here
    main(sys.argv[1:])
