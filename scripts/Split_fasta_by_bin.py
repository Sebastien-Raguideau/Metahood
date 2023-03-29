#!/usr/bin/env python
from __future__ import division
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import defaultdict
import argparse,os
import resource
import sys


def scheduler(method,List_Batchs,file,restrict,bin_composition,path_output,folder):
    for Batch in List_Batchs :
        if restrict :
            set_bins=set(Batch)&set(restrict)
        else :
            set_bins=set(Batch)
        method(file,bin_composition,set_bins,path_output,folder)


def split_fasta(fasta_file,bin_composition,set_bins,output,folder):
    # get bin definition 
    Dico_contigs_bin={line.rstrip().split(",")[0]:line.rstrip().split(",")[1] for line in open(bin_composition) if line.rstrip().split(",")[1] in set_bins}
    # open one file for writting, by bin
    Dico_bin_Handle={}
    ext=fasta_file.split(".")[-1]
    for bins in set(Dico_contigs_bin.values()) :
        if folder :
            output_folder=output+"/Bin_%s/"%bins
        else :
            output_folder=output
        os.system("mkdir -p %s"%output_folder)
        Dico_bin_Handle[bins]=open("%s/Bin_%s.%s"%(output_folder,bins,ext),"w")

    # go through all element of the the fasta file and write it to the rigth bin file
    for header,seq in SimpleFastaParser(open(fasta_file)) :
        contig1 = header.split()[0]
        contig2 = "_".join(contig1.split("_")[:-1])
        # deal with the fact that sometimes we get orfs and not contigs. 
        if contig1 in Dico_contigs_bin :
            Dico_bin_Handle[Dico_contigs_bin[contig1]].write(">%s\n%s\n"%(header,seq))
        if contig2 in Dico_contigs_bin :
            Dico_bin_Handle[Dico_contigs_bin[contig2]].write(">%s\n%s\n"%(header,seq))

    # close all the files
    for handle in Dico_bin_Handle.values() :
        handle.close()



def split_annotation(annotation_file,bin_composition,set_bins,output,folder):
    # get bin definition 
    Dico_contigs_bin={line.rstrip().split(",")[0]:line.rstrip().split(",")[1] for line in open(bin_composition) if line.rstrip().split(",")[1] in set_bins}
    # open one file for writting, by bin
    Dico_bin_Handle={}
    ext=annotation_file.split("/")[-1].replace("contigs_","")
    for bins in set(Dico_contigs_bin.values()) :
        if folder :
            output_folder=output+"/Bin_%s"%bins
        else :
            output_folder=output
        os.system("mkdir -p %s"%output_folder)
        Dico_bin_Handle[bins]=open("%s/Bin_%s_%s"%(output_folder,bins,ext),"w")

    # go through all element of the the fasta file and write it to the rigth bin file
    for line in open(annotation_file) :
        contig1 = line.rstrip().split("\t")[0]
        contig2 = "_".join(contig1.split("_")[:-1])
        for contig in [contig1,contig2] :
            if contig in Dico_contigs_bin :
                Dico_bin_Handle[Dico_contigs_bin[contig]].write(line)

    # close all the files
    for handle in Dico_bin_Handle.values() :
        handle.close()




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bin_composition", help="csv file giving bin composition, output of concoct, first column=contig id, second is bin number")
    parser.add_argument("path_output", help="specify the place you want all these bin folder to be put")
    parser.add_argument("--fasta",nargs='*', help="fasta files containing contigs or orfs binned")
    parser.add_argument("--annotation",nargs='*', help="annotation refering to sequences binned, must be .tsv, first column is the orf id, it should be possible to obtain contig id from removing last underscore suffix.")
    parser.add_argument("-l",help="file with one bin number by line : restrict to a list of specifics bins ")
    parser.add_argument("--folder",help="alternative output option : create one folder for each bin", action = 'store_true')
    args = parser.parse_args()

    # parse argument
    fasta_files=[]
    annotation_files=[]
    if args.fasta :
        fasta_files=[file for file in args.fasta]
    if args.annotation :
        annotation_files=[file for file in args.annotation]
    if annotation_files+fasta_files==[]:
        sys.exit('no fasta of annotation to split')
    bin_composition=args.bin_composition
    path_output=args.path_output
    folder=args.folder
    if args.l :
        restrict={line.rstrip() for line in open(args.l)}
    else :
        restrict=set()

    #---------- Batch strategy ----------------
    # there is a hard limit of number of handles which can be open
    max_nb_handles=resource.getrlimit(resource.RLIMIT_NOFILE)[0]
    max_nb_handles=max_nb_handles-100
    List_all_bins=list({line.rstrip().split(",")[1] for index,line in enumerate(open(bin_composition)) if index>0})
    Nb_bins=len(List_all_bins)
    nb_batch=Nb_bins//max_nb_handles+1
    List_Batchs=[List_all_bins[batch*max_nb_handles:min((batch+1)*max_nb_handles,Nb_bins)] for batch in range(nb_batch)]
    # split fasta files
    for file in fasta_files:
        scheduler(split_fasta,List_Batchs,file,restrict,bin_composition,path_output,folder)
    # split annotation files
    for file in annotation_files :
        scheduler(split_annotation,List_Batchs,file,restrict,bin_composition,path_output,folder)


