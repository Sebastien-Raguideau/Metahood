#!/usr/bin/env python3
# -*- coding: latin-1 -*-

from collections import defaultdict,Counter
from Bio.SeqIO.FastaIO import *
import numpy as np
import argparse

def get_linear_path(gfa_file) :
	def reverse_sign(sign):
		sign="+"*(sign=="-")+"-"*(sign=="+")
		return sign 
	Handle=open(gfa_file)
	Dico_connections=defaultdict(set)
	Dico_contig_line={}
	Dico_Path_lines=defaultdict(list)
	Set_Contigs=set([])
	List_edges=[]
	Self_loop_0=set([])
	Dico_contig_degree=defaultdict(lambda: [set([]),set([])])
	Dico_bicontigs_orientation=defaultdict()
	for Line in Handle : 
		Line2=Line.rstrip().split("\t")
		if Line2[0]=="S" : 
			# [S,name,sequence,LN:i:117,KC:i:85209,CL:z:#5cb9be]
			Dico_contig_line[Line2[1]]=Line
		if Line2[0]=="L" :
			# ['L','2','+','27','+','99M']
			(_,Contig1,Sign1,Contig2,Sign2,_)=Line2
			Dico_bicontigs_orientation[tuple(sorted([Contig1,Contig2]))]=Sign1==Sign2
			Dico_connections[(Contig1,Contig2)]|=set([(Sign1,Sign2)])
			Dico_connections[(Contig2,Contig1)]|=set([(reverse_sign(Sign2),reverse_sign(Sign1))])
			if Contig2==Contig1 :
				Self_loop_0|=set([Contig1])
			if Sign1 == "+" :
				Dico_contig_degree[Contig1][1]|={Contig2}
			else :
				Dico_contig_degree[Contig1][0]|={Contig2}
			if Sign2 == "+" :
				Dico_contig_degree[Contig2][0]|={Contig1}
			else :
				Dico_contig_degree[Contig2][1]|={Contig1}
	Handle.close()
	Dico_connections={key:values for key,values in Dico_connections.items() if len(values)>1}
	Self_loop_1=set(map(lambda x:tuple(sorted(x)),Dico_connections.keys()))
	return Dico_contig_degree,Dico_contig_line,Self_loop_0,Self_loop_1,Dico_bicontigs_orientation

def Crawl_through_subgraph_node(node_dir,Set_node,Dico_contig_degree) :
	(node,Direction)=node_dir
	List_New_node=Dico_contig_degree[node][Direction]
	if List_New_node :
		for New_node in List_New_node :
			Direction = node in Dico_contig_degree[New_node][0]
			if (New_node,Direction) not in Set_node : 
				Set_node.add((New_node,Direction))
				Crawl_through_subgraph_node((New_node,Direction),Set_node,Dico_contig_degree)


def write_subgraph(Set_nodes,Gfa_file,output_file) :
	Handle=open(output_file,"w")
	for Line in open(Gfa_file) : 
		if Line[0]=="S" : 
			# [S,name,sequence,LN:i:117,KC:i:85209,CL:z:#5cb9be]
			Contig=Line.rstrip().split("\t")[1]
			if Contig in Set_nodes :
				Handle.write(Line)
		if Line[0]=="L" :
		# ['L','2','+','27','+','99M']
			(_,Contig1,Sign1,Contig2,Sign2,_)=Line.rstrip().split("\t")
			if len({Contig1,Contig2}&Set_nodes)==2 :
					Handle.write(Line)
	Handle.close()


def contiguous_nodes_extraction(List_nodes_ini,Dico_contig_degree) :
	Set_contigs2=set()
	for contig in List_nodes_ini :
		for Direction in [0,1] :
			local_set={(contig,Direction)}
			Crawl_through_subgraph_node((contig,Direction),local_set,Dico_contig_degree)
			if len(local_set)>1 :
				Set_contigs2|=local_set
	return {orf for orf,direction in Set_contigs2}



def graph_selection(List_nodes,Gfa_file,ouptut_file) :
	# get graph infos
	Dico_contig_degree,Dico_contig_line,Self_loop_0,Self_loop_1,Dico_bicontigs_orientation=get_linear_path(Gfa_file)

	# check nodes from List_nodes are actualy in the graph
	Fraction=len(set(List_nodes)&set(Dico_contig_line.keys()))/float(len(set(List_nodes)))
	if Fraction!=1 :
		"Warning %2.2f %% of the contigs were found in the graph, check the contigs inputed (for instance check they are not ORFs) "% (Fraction)

	# select nodes contiguous to the ones in List_nodes
	Set_nodes_selected=contiguous_nodes_extraction(List_nodes,Dico_contig_degree)

	# read the initial graph and only retain line concerning Set_nodes_selected for writing
	write_subgraph(Set_nodes_selected,Gfa_file,ouptut_file)


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
