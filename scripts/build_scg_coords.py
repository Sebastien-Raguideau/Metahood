#!/usr/bin/env python3

from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", help="SCG bed")
    parser.add_argument("-d", help="bin definition")
    parser.add_argument("-f", help="SCG file")
    parser.add_argument("-o", help="output")
    parser.add_argument("-w", help="wildcards.bin")
    args = parser.parse_args()
    
    scg_bed=args.b
    bin_def=args.d
    scg_file=args.f
    out_coor=args.o
    wild_bin=args.w

    # contigs specific to the bin
    set_contigs={line.rstrip().split(",")[0] for line in open(bin_def) if (line.rstrip().split(",")[1]==wild_bin.replace("Bin_",""))}
    # NODE_11_length_160076_cov_189.092050_50 COG0552 strand=+
    # Cog annotation of orfs in header of the file
    dict_scg = {header.split(" ")[0]: header.split(" ") for header, seq in sfp(open(scg_file))}
    with open(out_coor, "w") as handle :
        for line in open(scg_bed):
            contig, start, end, orf = line.rstrip().split("\t")
            if contig in set_contigs:
                cog, strand = dict_scg[orf][1:]
                cog, strand = dict_scg[orf][1:]
                strand = str((strand == "strand=+")*1-(strand == "strand=-")*1)
                contig = "_".join(orf.split("_")[:-1])
                handle.write(",".join([cog, contig, start, end, orf, strand])+"\n")

