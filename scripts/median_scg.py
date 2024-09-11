#!/usr/bin/env python3

import numpy as np
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", help="input cogs tsv")
    parser.add_argument("-n", help="input nuc tsv")
    parser.add_argument("-o", help="output tsv")
    parser.add_argument("-s", help="SCG_DATA")
    args = parser.parse_args()
    
    cov=args.c
    nuc=args.n
    out_file=args.o
    SCG_DATA=args.s


    set_SCG = {cog.rstrip() for cog in open(SCG_DATA+"/scg_cogs_min0.97_max1.03_unique_genera.txt")}
    List_profile = []
    # compute median of SCG
    with open(cov) as handle:
        samples = next(handle).rstrip().split("\t")[1:]
        for index,line in enumerate(handle):
            split_line = line.rstrip().split("\t")
            cog = split_line[0]
            if cog in set_SCG:
                List_profile.append([float(element) for element in split_line[1:]])

    if not List_profile:
        scg_norm = np.array([0 for s in samples])
    else:
        scg_norm=np.median(List_profile, axis=0)

    # get previous normalisation
    with open(nuc) as handle:
        samples_local=next(handle).rstrip().split("\t")[1:]
        nuc_norm=next(handle).rstrip().split("\t")[1:]

    with open(out_file,"w") as handle:
        handle.write("Normalisation\t"+"\t".join(samples_local)+"\n")
        handle.write("Nucleotides\t"+"\t".join(nuc_norm)+"\n")
        handle.write("median_scg\t"+"\t".join([str(scg_norm[samples.index(sample)]) for sample in samples_local])+"\n")

