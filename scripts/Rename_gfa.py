#!/usr/bin/env python3
# -*- coding: latin-1 -*-
import argparse
import os
from subprocess import Popen,PIPE
import random 
from Bio.SeqIO.FastaIO import *


def Rewrite_gfa(fasta_file,gfa_file,output) :
	# get contig name since megahit_toolkit delete them (grr)
	def Get_contig_id(Handle):
		Ordered_Contig_id=[]
		for title, sequence in SimpleFastaParser(Handle) :
			Ordered_Contig_id.append(title.split(' ')[0])
		return Ordered_Contig_id
	Handle=open(fasta_file)
	Dico_Order_Contig={str(Order+1):Name for Order,Name in enumerate(Get_contig_id(Handle))}
	Handle.close()
	Handle=open(gfa_file)
	NewGfa=""
	for line in Handle :
		line=line.rstrip().split('\t')
		if line[0]=="S" :
			line[-1]+='\n'
			line[1]=Dico_Order_Contig[line[1]]
		elif line[0]=="L":
			line[1]=Dico_Order_Contig[line[1]]
			line[3]=Dico_Order_Contig[line[3]]
			line.append('\n')
		NewGfa+="\t".join(line)
	Handle.close()
	H=open(output,'w')
	H.write(NewGfa)
	H.close()
	return Dico_Order_Contig.values()

def test_random_contig(Ordered_Contig_id,fasta_file,gfa_file,output) :
	Contig=random.choice(list(Ordered_Contig_id))
	# get seq in intial fasta file
	process = Popen(['grep -A1 "'+Contig+'" '+fasta_file], stdout=PIPE, stderr=PIPE,shell=True,executable='/bin/bash')
	Seq_fasta=process.communicate()[0].decode('utf-8').split("\n")[1]
	# get seq in translated gfa file 
	process = Popen(["grep $'"+Contig+"\t' " +output], stdout=PIPE, stderr=PIPE,shell=True,executable='/bin/bash')
	Seq_gfa=process.communicate()[0].decode('utf-8').split("\t")[2]
	if Seq_fasta!=Seq_gfa :
		print("warning renaming is wrong, "+Contig +"is different in both files")

	
def main(gfa_file,fasta_file,output) :
	Ordered_Contig_id=Rewrite_gfa(fasta_file,gfa_file,output)
	test_random_contig(Ordered_Contig_id,fasta_file,gfa_file,output)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("gfa_file", help="gfa file you want to rename") 
	parser.add_argument("fasta_file", help="initial fasta file containing the right names, contigs should be ordered in the same way")
	parser.add_argument("output", help="output filename")
	args = parser.parse_args()
	gfa_file=args.gfa_file
	fasta_file=args.fasta_file
	output=args.output
	main(gfa_file,fasta_file,output)
