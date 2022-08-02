#!/usr/bin/env python3
import argparse
from collections import defaultdict,Counter 
import numpy as np

def main(table,cog_hmm,out_file):
    # get cog definition in terms of pfam/TIGR
    hmm_to_cog = {el:line.rstrip().split("\t")[0] for line in open(cog_hmm) if line.rstrip().split("\t")[2]=="fine" for el in line.rstrip().split("\t")[3:]}


    # get best hits hmm
    orf_hmm = defaultdict(list)
    for line in open(table):
        if line[0]!="#":
            orf_hmm[line.rstrip().split()[0]].append([line.rstrip().split()[4],line.rstrip().split()[7]])

    best_hit_hmm = {orf:max(hmm,key=lambda x:float(x[1]))[0] for orf,hmm in orf_hmm.items()}
    best_hit_hmm = {orf:hmm.split(".")[0] for orf,hmm in best_hit_hmm.items()}

    # get orf annotation
    orf_to_cog = {orf:hmm_to_cog[hmm] for orf,hmm in best_hit_hmm.items() if hmm in hmm_to_cog}

    # output:
    with open(out_file,"w") as handle:
        handle.write("Querry\tSubject\n")
        handle.writelines("%s\t%s\n"%(orf,cog) for orf,cog in orf_to_cog.items())





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("table",help='hmm search output table')    
    parser.add_argument('--cog_hmm',help = ('mapping of cog to hmm profiles'))
    parser.add_argument('output')    
    args = parser.parse_args()
    main(args.table,args.cog_hmm,args.output)
