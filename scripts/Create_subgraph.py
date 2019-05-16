#!/usr/bin/env python3
# -*- coding: latin-1 -*-
import argparse
from collections import defaultdict

def Get_graph_info(Gfa_file) :
	Dico_contig_links=defaultdict(set)
	Dico_contig_neighbors=defaultdict(set)
	for Line in open(Gfa_file) : 
		if Line[0]=="L" :
			# ['L','2','+','27','+','99M']
			(_,Contig1,Sign1,Contig2,Sign2,_)=Line.rstrip().split("\t")
			Set_contigs=set([Contig1,Contig2])
			Dico_contig_links[Contig1]|=set([Line])
			Dico_contig_links[Contig2]|=set([Line])
			Dico_contig_neighbors[Contig1]|=set([Contig2])
			Dico_contig_neighbors[Contig2]|=set([Contig1])
	return Dico_contig_links,Dico_contig_neighbors

def Get_seq_info(Gfa_file,Set_contigs) :
	Dico_contig_line={}
	for Line in open(Gfa_file) : 
		if Line[0]=="S" : 
			Contig=Line.rstrip().split("\t")[1]
			if Contig in Set_contigs :
				Dico_contig_line[Line.rstrip().split("\t")[1]]=Line
			# [S,name,sequence,LN:i:117,KC:i:85209,CL:z:#5cb9be]
	return Dico_contig_line

def Write_gfa(Dico_contig_line,Dico_contig_links,ouptut_file):
	Edges="".join({edge for contig in Dico_contig_line for edge in Dico_contig_links[contig]})
	Handle=open(ouptut_file,"w")
	Handle.write("".join(Dico_contig_line.values()))
	Handle.write(Edges)
	Handle.close()

def graph_selection(List_nodes,Gfa_file,ouptut_file) :
	Dico_contig_links,Dico_contig_neighbors=Get_graph_info(Gfa_file)
	New_contig_set=set([])
	To_Process=set(List_nodes)
	Fraction=len(To_Process&set(Dico_contig_links.keys()))/float(len(To_Process))
	if Fraction!=1 :
		"Warning %2.2f %% contigs were found to be conencted in the initial graph, check the contigs inputed (for instance check they are not ORFs) "% (Fraction)
	while To_Process :
		New_contig_set|=To_Process
		To_Process|={contig for contig_to_process in To_Process for contig in Dico_contig_neighbors[contig_to_process]}
		To_Process-=New_contig_set
	Dico_contig_line=Get_seq_info(Gfa_file,New_contig_set)
	Write_gfa(Dico_contig_line,Dico_contig_links,ouptut_file)



if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("Annotation_file", help="Annotation file, first column are ORFs, corresponding contig can be obtained by just stripping the last underscore")
	parser.add_argument("Gfa_file", help="Gfa file you want to extract a subgraph from")
	parser.add_argument("-o", help="output file name")
	args = parser.parse_args()
	Annotation_file=args.Annotation_file
	List_nodes=list({"_".join(orf.split("\t")[0].split('_')[:-1]) for orf in open(Annotation_file)})
	Gfa_file=args.Gfa_file
	if args.o :
		ouptut_file=args.o
	else :
		ouptut_file=Gfa_file.replace(".gfa","_"+Node_filename+".gfa")
	graph_selection(List_nodes,Gfa_file,ouptut_file)

# Node_filename="List_integron_loose.tsv"
# Gfa_file="./Assembly_final.gfa"
# ouptut_file="test"
# List_nodes=[contig.split(".")[0] for contig in open(Node_filename).read().rstrip().split("\t")]
# Missing_contigs=New_contig_set-set(Dico_contig_line.keys())
