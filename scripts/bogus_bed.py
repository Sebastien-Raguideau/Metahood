#!/usr/bin/env python3

from Bio.SeqIO.FastaIO import SimpleFastaParser as SFP
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input contig")
    parser.add_argument("-o", help="output bed")
    args = parser.parse_args()
    input_contig=args.i
    output_bed=args.o
 
    handle=open(output_bed,"w")
    for header,seq in SFP(open(input_contig)) :
        name=header.split(" ")[0]
        handle.write("\t".join([name,"0",str(len(seq)),name+"\n"]))
    handle.close()
