#!/usr/bin/env python3

from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", help="orfs bed")
    parser.add_argument("-f", help="SCG file")
    parser.add_argument("-o", help="output SCG bed")
    args = parser.parse_args()
    orfs_bed=args.b
    scg_file=args.f
    scg_bedfile=args.o


    seq_names = {name.split()[0] for name, seq in sfp(open(scg_file))}
    with open(scg_bedfile, 'w') as scg_bed:
        for line in open(orfs_bed):
            if line.split()[3] in seq_names:
               scg_bed.write(line)

