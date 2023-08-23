#!/usr/bin/env python
import os 
from os.path import realpath,abspath,basename
from Bio.SeqIO.FastaIO import *
from multiprocessing import Pool
from subprocess import Popen, PIPE
from collections import defaultdict
import argparse
import time 

def split_fasta(fasta_file,Temp_location,nb_chunks):
    Sorted_Names=[]
    Dico_genome_seq={}
    Dico_genome_len={}
    for header,seq in SimpleFastaParser(open(fasta_file)) :
        Sorted_Names.append(header)
        Dico_genome_seq[header]=seq
        Dico_genome_len[header]=len(seq)
    Total_length = sum(Dico_genome_len.values())
    Chunk_size = Total_length/float(nb_chunks)
    os.system("mkdir -p "+Temp_location)

    # schedule using mean size reassessed as we define batchsize
    fasta_path = "%s/Batch"%Temp_location
    Current_filename=lambda x:"%s_%s"%(fasta_path,str(x))
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

def prodigal(file):
    process = Popen(["prodigal -i %s -a %s.faa -d %s.fna -f gff -p meta -o %s.gff > %s_prodigal.out 2>&1"%(file,file,file,file,file)], stdout=PIPE, stderr=PIPE,shell=True)
    return process.communicate()


def main(fasta_file,Temp_location,threads,nb_chunks,output_path) :

    # split fasta file in temp_split/Batch
    split_fasta(fasta_file,Temp_location,nb_chunks)

    # list files
    files = [realpath("%s/Batch_%s"%(Temp_location,nb)) for nb in range(nb_chunks)]

    # parrallelize stuf
    pool = Pool(threads)
    result = pool.map(prodigal,files)

    # gather batch in unique files 
    new_name  = "%s/%s"%(output_path,".".join(basename(fasta_file).split('.')[:-1]))
    with open(("%s.faa"%new_name),"w") as h_faa, open(("%s.fna"%new_name),"w") as h_fna, open(("%s.gff"%new_name),"w") as h_gff:
        h_gff.write(open("%s/Batch_0.gff"%Temp_location).readline())
        for index in range(nb_chunks) :
            h_faa.write(open("%s/Batch_%s.faa"%(Temp_location,index)).read())
            h_fna.write(open("%s/Batch_%s.fna"%(Temp_location,index)).read())
            tmp=open("%s/Batch_%s.gff"%(Temp_location,index))
            _=tmp.readline()
            h_gff.write(tmp.read())

    time.sleep(10)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("n", help="Number of prodigal tasks")
    parser.add_argument("-s", help="Number of chunks, default is same as number of tasks")
    parser.add_argument("-T", help="Temp file location ",default="./split_annotation")
    parser.add_argument("f", help="fasta_file")
    parser.add_argument("-o", help="output Directory",default='./')
    args = parser.parse_args()
    fasta_file=args.f
    threads=int(args.n)
    Temp_location=args.T
    Temp_location=Temp_location+(not Temp_location[-1]=="/")*"/"
    if args.s :
        nb_chunks=int(args.s)
    else:
        nb_chunks=threads
    output_path=args.o
    main(fasta_file,Temp_location,threads,nb_chunks,output_path)




