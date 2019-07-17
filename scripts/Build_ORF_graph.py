#!/usr/bin/env python3
# -*- coding: latin-1 -*-
import argparse
from collections import defaultdict
import re


def Get_contig_description(bed_file,Dico_contig_lines) :
	Dico_Contigs_ORFs=defaultdict(list)
	Dico_ORF_Seq={}
	for line in open(bed_file) :
		Contig,start,end,ORF=line.rstrip().split("\t")
		if Contig in Dico_contig_lines :
			Seq=Dico_contig_lines[Contig].split("\t")[2]
			Dico_ORF_Seq[ORF]=Seq[int(start):int(end)]
			Dico_Contigs_ORFs[Contig].append(ORF)
	for Contig in Dico_contig_lines :
		if Contig not in Dico_Contigs_ORFs :
			Dico_Contigs_ORFs[Contig]=[Contig]
			Seq=Dico_contig_lines[Contig].split("\t")[2]
			Dico_ORF_Seq[Contig]=Seq
	return Dico_Contigs_ORFs,Dico_ORF_Seq

def Parse_GFA(Gfa_file) :
	Dico_contig_line={}
	List_edges_ini=[]
	for Line in open(Gfa_file) : 
		if Line[0]=="S" : 
			Dico_contig_line[Line.rstrip().split("\t")[1]]=Line
		if Line[0]=="L" :
			List_edges_ini.append(Line)
	return Dico_contig_line,List_edges_ini

def Write_gfa_ORF(output,Dico_ORF_Seq,Dico_contig_lines,Dico_Contigs_ORFs,List_edges_ini) : 
	List_match=list({line.split('\t')[5][:-1] for line in List_edges_ini})
	if len(List_match)>1 :
		print("more than one Kmer length : "+"\t".join(List_match)) 
	else :
		kmer_len=int(List_match[0])
	List_lines=[]
	List_edges=[]
	for Contig,Line in Dico_contig_lines.items() :
		# [S,name,sequence,LN:i:117,KC:i:85209]
		KC_contig=0
		Contig_len=len(Line.split("\t")[2])
		Match_KC=re.search("KC:i:(\d*)",Line).groups()
		if Match_KC :
			KC_contig=int(Match_KC[0])
		for index,ORF in enumerate(Dico_Contigs_ORFs[Contig]) :
			Seq=Dico_ORF_Seq[ORF]
			ORF_Len=len(Seq)
			LN="LN:i:"+str(ORF_Len)
			KC="KC:i:"+ str((KC_contig/(Contig_len-kmer_len))*ORF_Len)
			List_lines.append("\t".join(["S",ORF,str(Seq),LN,KC]))
			if index!=0 :
				List_edges.append("\t".join(["L",Previous_ORF,"+",ORF,"+","1M"]))
			Previous_ORF=ORF
	for Line in List_edges_ini :
		(_,Contig1,sgn1,Contig2,sgn2,M)=Line.rstrip().split("\t")
		if sgn1=="+" :
			ORF1=Dico_Contigs_ORFs[Contig1][-1]
		else :
			ORF1=Dico_Contigs_ORFs[Contig1][0]
		if sgn2=="+" :
			ORF2=Dico_Contigs_ORFs[Contig2][0]
		else :
			ORF2=Dico_Contigs_ORFs[Contig2][-1]
		List_edges.append("\t".join(["L",ORF1,sgn1,ORF2,sgn2,M]))
	Handle=open(output,"w")
	Handle.write("\n".join(List_lines))
	Handle.write("\n".join(List_edges))
	Handle.close()

def main(bed_file,Gfa_file,output) :
	Dico_contig_lines,List_edges_ini=Parse_GFA(Gfa_file)
	Dico_Contigs_ORFs,Dico_ORF_Seq=Get_contig_description(bed_file,Dico_contig_lines)
	Write_gfa_ORF(output,Dico_ORF_Seq,Dico_contig_lines,Dico_Contigs_ORFs,List_edges_ini)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("bed_file", help="tsv file wiht field contig\tstart\tend\tORF ")
	parser.add_argument("Gfa_file", help="GFA file you want to convert")
	parser.add_argument("output", help="output file")
	args = parser.parse_args()
	bed_file=args.bed_file
	Gfa_file=args.Gfa_file
	output=args.output
	main(bed_file,Gfa_file,output)






# Gfa_file="AMR_reduced2.gfa"
# bed2_file="Annotation/Contigs_F500_C10K.bed2"








#################### Stupid stuff #####################


# def Parse_gfa_light(gfa_file,G):
# 	def reverse_sign(sign):
# 		sign="+"*(sign=="-")+"-"*(sign=="+")
# 		return sign 
# 	def reverse_line(Line) :
# 		Line2=[Line[0],Line[3],reverse_sign(Line[4]),Line[1],reverse_sign(Line[2]),Line[5]]
# 		return Line2
# 	Dico_node_infos={}
# 	Handle=open(gfa_file)
# 	for Line in Handle : 
# 		Line=Line.rstrip().split("\t")
# 		if Line[0]=="S" : 
# 			# [S,name,sequence,LN:i:117,KC:i:85209,CL:z:#5cb9be]
# 			(name,seq_plus)=(Line[1],Seq(Line[2]))
# 			G.add_node(name+"+",Name=name+"+",Cov_len="\t".join(Line[3:]))
# 			G.add_node(name+"-",Name=name+"-",Cov_len="\t".join(Line[3:]))
# 		if Line[0]=="L" :
# 			# ['L','2','+','27','+','99M']
# 			G.add_edge("".join(Line[1:3]),"".join(Line[3:5]))
# 			G["".join(Line[1:3])]["".join(Line[3:5])]['Struct']="\t".join(Line)
# 			G.add_edge("".join([Line[3],reverse_sign(Line[4])]),"".join([Line[1],reverse_sign(Line[2])]))
# 			G[Line[3]+reverse_sign(Line[4])][Line[1]+reverse_sign(Line[2])]['Struct']="\t".join(reverse_line(Line))
# 	Handle.close()

# def Get_back_Additionals_contigs(gfa_file,set_missing_contigs,G):
# 	Missing=0
# 	Handle=open(gfa_file)
# 	for Line in Handle : 
# 		Line=Line.rstrip().split("\t")
# 		if Line[0]=="S" : 
# 			# [S,name,sequence,LN:i:117,KC:i:85209,CL:z:#5cb9be]
# 			(name,seq_plus)=(Line[1],Seq(Line[2]))
# 			if name in set_missing_contigs :
# 				G.add_node(name+"+",Name=name+"+",Cov_len="\t".join(Line[3:]),Seq=seq_plus)
# 				G.add_node(name+"-",Name=name+"-",Cov_len="\t".join(Line[3:]))
# 				Missing+=1
# 				if Missing==len(set_missing_contigs) :
# 					break
# 		if Line[0]=="L" :
# 			break
# 	Handle.close()



# def Get_G_reduced(Gfa_file,List_Contigs,Gfa_source) :
# 	G=nx.Graph()
# 	Parse_gfa_light(Gfa_file,G)
# 	if List_Contigs :
# 		List_Contigs=[contig.split(".")[0] for contig in  List_Contigs]
# 		if Gfa_source :
# 			Set_Missing_contigs=set(List_Contigs)-set([node[:-1] for node in G.nodes()])
# 			if Set_Missing_contigs :
# 				Get_back_Additionals_contigs(Gfa_source,Set_Missing_contigs,G)
# 		Dico_nodes_subgraphs={node:subgraph for subgraph in list(nx.connected_component_subgraphs(G)) for node in subgraph.nodes()}
# 		List_nodes=[Contigs+orientation for Contigs in List_Contigs for orientation in ["+","-"] if Contigs+orientation in Dico_nodes_subgraphs]
# 		List_reduced_subgraphs=list({Dico_nodes_subgraphs[node] for node in List_nodes})
# 		G_reduced=List_reduced_subgraphs[0]
# 		for subgraph in List_reduced_subgraphs[1:] :
# 			G_reduced=nx.compose(G_reduced, subgraph)
# 		return G_reduced
# 	else :
# 		return G/

# def GFA_get_node_seq(gfa_file,List_contig):
# 	Set_Contigs=set(List_contig)
# 	Dico_contig_seq={}
# 	Handle=open(gfa_file)
# 	for Line in Handle : 
# 		Line=Line.rstrip().split("\t")
# 		if Line[0]=="S" : 
# 			# [S,name,sequence,LN:i:117,KC:i:85209,CL:z:#5cb9be]
# 			(name,seq_plus)=(Line[1],Seq(Line[2]))
# 			if name in Set_Contigs :
# 				Dico_contig_seq[name]=seq_plus
# 	return Dico_contig_seq


# def get_List_contigs_from_CARD(Card_annotation) :
# 	#get corresponding contigs
# 	Dico_contig_CARD=defaultdict(list)
# 	Handle=open(Card_annotation)
# 	Header=Handle.next()
# 	List_ORFs=[]
# 	for line in Handle :
# 		AMR=line.split("\t")[1]
# 		ORF=line.split("\t")[0]
# 		List_ORFs.append(ORF)
# 		Contig="_".join(ORF.split("_")[:2])
# 		Dico_contig_CARD[Contig].append(AMR)
# 	Handle.close()
# 	Dico_contig_CARD_multiple={contig:list_amr for contig,list_amr in Dico_contig_CARD.items() if len(list_amr)>1}
# 	## get List of contig of interest : 
# 	List_Contigs=list(Dico_contig_CARD.keys())
# 	return List_Contigs



