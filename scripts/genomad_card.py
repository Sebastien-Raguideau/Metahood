#!/usr/bin/env python
import argparse
import numpy as np
from collections import defaultdict
from os.path import dirname


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("card", help="type of annotation")
    parser.add_argument("genomad", help="list of assemblies")
    parser.add_argument("output", help="output")
    args = parser.parse_args()
    card_file = args.card
    genomad_dir =  dirname(args.genomad)
    output = args.output
        
    card_contigs = {"_".join(line.rstrip().split("\t")[0].split("_")[:-1]):[line.rstrip().split("\t")[1],line.rstrip().split("\t")[0]] for index,line in enumerate(open(card_file))  if index>0}
    vir_contigs = {line.rstrip().split("\t")[0] for index,line in enumerate(open(f"{genomad_dir}/contigs_summary/contigs_virus_summary.tsv")) if index>0}
    plas_contigs = {line.rstrip().split("\t")[0] for index,line in enumerate(open(f"{genomad_dir}/contigs_summary/contigs_plasmid_summary.tsv")) if index>0}

    card_map = defaultdict(list)
    for cont,(aro,orf) in card_contigs.items():
        if cont in plas_contigs:
            card_map["%s|plasmid"%aro].append(orf)
        elif cont in vir_contigs:
            card_map["%s|virus"%aro].append(orf)
        else:
            card_map[aro].append(orf)

    with open(output,"w") as handle:
        for aro,orfs in card_map.items():
            handle.writelines("%s\t%s\n"%(aro,"\t".join(orfs)))







