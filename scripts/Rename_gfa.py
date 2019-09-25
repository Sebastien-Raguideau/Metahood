#!/usr/bin/env python
# -*- coding: latin-1 -*-
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
import argparse
import os

def Rewrite_gfa(fasta_file,gfa_file,output) :
	# get contig name since megahit_toolkit delete them (grr)
	# may look ineficient but takes only 10 min in total on 16G assembly
	with open(fasta_file) as Handle :
		Dico_Seq_name={seq:name.split(' ')[0] for name,seq in sfp(Handle)}

	NewGfa=[]
	with open(gfa_file) as Handle :
		Dict_old_to_new={}
		for line in Handle :
			line=line.rstrip().split('\t')
			if line[0]=="S" :
				Dict_old_to_new[line[1]]=Dico_Seq_name[line[2]]
				line[1]=Dico_Seq_name[line[2]]
				NewGfa.append("\t".join(line))

	with open(gfa_file) as Handle :
		for line in Handle :
			line=line.rstrip().split('\t')
			if line[0]=="L":
				line[1]=Dict_old_to_new[line[1]]
				line[3]=Dict_old_to_new[line[3]]
				NewGfa.append("\t".join(line))

	with open(output,'w') as H :
		H.write("\n".join(NewGfa))

def main(gfa_file,fasta_file,output) :
	Ordered_Contig_id=Rewrite_gfa(fasta_file,gfa_file,output)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("gfa_file", help="gfa file you want to rename") 
	parser.add_argument("fasta_file", help="initial fasta file containing the right names")
	parser.add_argument("output", help="output filename")
	args = parser.parse_args()
	gfa_file=args.gfa_file
	fasta_file=args.fasta_file
	output=args.output
	main(gfa_file,fasta_file,output)
