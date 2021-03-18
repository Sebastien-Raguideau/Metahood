#!/usr/bin/env python3

from Bio.SeqIO.FastaIO import SimpleFastaParser as SFP
import argparse
import glob

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", help="Path to write to")
    parser.add_argument("-w", help="Wildcard group")
    args = parser.parse_args()
    out_file=args.o
    wild_grp=args.w

    List_bins=glob.glob(wild_grp+"/binning/metabat2/bins/Bin*.fa")
    Handle=open(out_file,"w")
    Handle.write("contig_id,0\n")
    for file in List_bins :
        bin_name=file.split("Bin_")[-1].split('.fa')[0]
        for name,seq in SFP(open(file)) :
            Handle.write(",".join([name,bin_name])+"\n")
    Handle.close()
