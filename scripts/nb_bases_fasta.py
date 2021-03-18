#!/usr/bin/env python3

from Bio.SeqIO.FastaIO import SimpleFastaParser as SFP
import argparse
import gzip

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input")
    parser.add_argument("-o", help="output")
    args = parser.parse_args()
    input=args.i
    output=args.o
    
    dict_sample_nb=defaultdict(list)
    for file in input:
        sample=basename(dirname(file))
        if file.endswith(".gz") :
            handle = gzip.open(file,"rt")
        else :
            handle = open(file)
        nb = sum([len(seq) for header,seq in sfp(handle)])
        dict_sample_nb[sample].append(nb)
    dict_sample_nb={sample:str(sum(list_nb)) for sample,list_nb in dict_sample_nb.items()}
    with open(output[0],"w") as handle:
        handle.write("Normalisation\t"+"\t".join(dict_sample_nb.keys())+"\n")
        handle.write("Nb_nucleotides\t"+"\t".join(dict_sample_nb.values())+"\n")
    handle.close()
