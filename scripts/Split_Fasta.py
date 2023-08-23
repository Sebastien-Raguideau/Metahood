#!/usr/bin/env python
# -*- coding: latin-1 -*-
from Bio.SeqIO.FastaIO import *
from collections import defaultdict,Counter
import numpy as np
import argparse
import os 

def split_fasta(fasta_file,nb_chunks,output) :
    Sorted_Names=[]
    Dico_genome_seq={}
    Dico_genome_len={}
    for header,seq in SimpleFastaParser(open(fasta_file)) :
        Sorted_Names.append(header)
        Dico_genome_seq[header]=seq
        Dico_genome_len[header]=len(seq)
    Total_length=sum(Dico_genome_len.values())
    Chunk_size=Total_length/float(nb_chunks)
    os.system("mkdir -p "+output)

    # schedule using mean size reassessed as we define batchsize
    extension="."+fasta_file.split(".")[-1]
    fasta_path=output+"/Batch"
    Current_filename=lambda x:fasta_path+"_"+str(x)+extension
    Temp_length=0
    num=0
    fname = Current_filename(num)
    batch_headers = defaultdict(list)
    for index,header in enumerate(Sorted_Names):
        if Temp_length>Chunk_size :
            Temp_length = 0
            num+=1
            fname = Current_filename(num)
            Chunk_size = Total_length/float(nb_chunks-num)

        # small file catchup for at least 1 seq per file
        if len(Dico_genome_len)-index<=(nb_chunks-num):
            for header in Sorted_Names[index:]:
                batch_headers[fname].append(header)
                num+=1
                fname = Current_filename(num)
            break

        # used for adaptative Chunk_size
        Len = Dico_genome_len[header]
        Temp_length+=Len
        Total_length-=Len

        # schedule which header in which file
        batch_headers[fname].append(header)

    # handle case it doesn't work for any reason
    if len(Dico_genome_seq)>nb_chunks:
        assert(num==nb_chunks),"splitting fasta failed, num=%s instead of %s, there is a bug"%(num,nb_chunks)

    # create the files
    for file,headers in batch_headers.items():
        with open(file,"w") as handle:
            handle.writelines(">%s\n%s\n"%(header,Dico_genome_seq[header]) for header in headers)

    for n in range(num+1,nb_chunks):
        os.system("touch %s"%Current_filename(n))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("Fasta_file", help="fasta file you want to split in chuncks")
    parser.add_argument("nb_chunks", help="nb of Chunks")
    parser.add_argument("-o", help="Name of the folder where you want to store the chunks",default="Split_")
    args = parser.parse_args()
    fasta_file=args.Fasta_file
    n=int(args.nb_chunks)
    output=args.o
    if output[-1]=="/" :
        output=output[:-1]
    if output=="Split_" :
        output+=".".join(fasta_file.split(".")[:-1])
    split_fasta(fasta_file,n,output)
