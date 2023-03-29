#!/usr/bin/env python
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
    for Line in Handle:
        Line2=Line.rstrip().split("\t")
        if Line2[0]=="S":
            # [S,name,sequence,LN:i:117,KC:i:85209,CL:z:#5cb9be]
            Dico_contig_line[Line2[1]]=Line
        if Line2[0]=="L":
            # ['L','2','+','27','+','99M']
            (_,Contig1,Sign1,Contig2,Sign2,_)=Line2
            Dico_bicontigs_orientation[tuple(sorted([Contig1,Contig2]))]=Sign1==Sign2
            Dico_connections[(Contig1,Contig2)]|=set([(Sign1,Sign2)])
            Dico_connections[(Contig2,Contig1)]|=set([(reverse_sign(Sign2),reverse_sign(Sign1))])
            if Contig2==Contig1 :
                Self_loop_0|=set([Contig1])
            if Sign1 == "+":
                Dico_contig_degree[Contig1][1]|={Contig2}
            else:
                Dico_contig_degree[Contig1][0]|={Contig2}
            if Sign2 == "+":
                Dico_contig_degree[Contig2][0]|={Contig1}
            else:
                Dico_contig_degree[Contig2][1]|={Contig1}
    Handle.close()
    Dico_connections={key:values for key,values in Dico_connections.items() if len(values)>1}
    Self_loop_1=set(map(lambda x:tuple(sorted(x)),Dico_connections.keys()))
    return Dico_contig_degree,Dico_contig_line,Self_loop_0,Self_loop_1,Dico_bicontigs_orientation

def Crawl_through_subgraph_node(node_Direction,Set_node,Dico_contig_degree,contig_len,dist,dist_max) :
    node,Direction=node_Direction
    List_New_node=Dico_contig_degree[node][Direction]
    if dist>dist_max:
    	return
    if List_New_node:
        for New_node in List_New_node:
            Direction = node in Dico_contig_degree[New_node][0]
            if (New_node,Direction) not in Set_node : 
                Set_node.add((New_node,Direction))
                dist+=contig_len[New_node]
                Crawl_through_subgraph_node((New_node,Direction),Set_node,Dico_contig_degree,contig_len,dist,dist_max)


def write_subgraph(Set_nodes,Gfa_file,output_file) :
    Handle=open(output_file,"w")
    for Line in open(Gfa_file):
        if Line[0]=="S":
            # [S,name,sequence,LN:i:117,KC:i:85209,CL:z:#5cb9be]
            Contig=Line.rstrip().split("\t")[1]
            if Contig in Set_nodes :
                Handle.write(Line)
        if Line[0]=="L":
        # ['L','2','+','27','+','99M']
            (_,Contig1,Sign1,Contig2,Sign2,_)=Line.rstrip().split("\t")
            if len({Contig1,Contig2}&Set_nodes)==2:
                    Handle.write(Line)
    Handle.close()


def contiguous_nodes_extraction(List_nodes_ini,Dico_contig_degree,contig_len,dist_max) :
    Set_contigs2=set()
    for contig in List_nodes_ini :
        for Direction in [0,1] :
            local_set={(contig,Direction)}
            Crawl_through_subgraph_node((contig,Direction),local_set,Dico_contig_degree,contig_len,0,dist_max)
            if len(local_set)>1 :
                Set_contigs2|=local_set
    return {orf for orf,direction in Set_contigs2}



def graph_selection(List_nodes,Gfa_file,mode,dist,ouptut_file):
    # get graph infos
    Dico_contig_degree,Dico_contig_line,Self_loop_0,Self_loop_1,Dico_bicontigs_orientation=get_linear_path(Gfa_file)

    # check nodes from List_nodes are actualy in the graph
    Fraction=len(set(List_nodes)&set(Dico_contig_line.keys()))/float(len(set(List_nodes)))
    if Fraction!=1 :
        "Warning %2.2f %% of the contigs were found in the graph, check the contigs inputed (for instance check they are not ORFs) "% (Fraction)

    if mode =="contiguous" :
        # select nodes contiguous to the ones in List_nodes
        contig_len = {contig:len(line.rstrip().split("\t")[2]) for contig,line in Dico_contig_line.items()}
        Set_nodes_selected = contiguous_nodes_extraction(List_nodes,Dico_contig_degree,contig_len,dist)

    if mode=="component" :
        # select nodes in the same components to the ones in List_nodes
        To_Process=set(List_nodes)
        # create contig to neighbors dict 
        Dico_contig_neighbors=defaultdict(set)
        for contig,(F,R) in Dico_contig_degree.items() :
            Dico_contig_neighbors[contig]=F|R
        # select nodeswith a breadth first approach 
        Set_nodes_selected=set()
        while To_Process :
            Set_nodes_selected|=To_Process
            To_Process|={contig for contig_to_process in To_Process for contig in Dico_contig_neighbors[contig_to_process]}
            To_Process-=Set_nodes_selected

    # read the initial graph and only retain line concerning Set_nodes_selected for writing
    Set_nodes_selected|=set(List_nodes)
    write_subgraph(Set_nodes_selected,Gfa_file,ouptut_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("Annotation_file", help="Annotation file, first column are ORFs, corresponding contig can be obtained by just stripping the last underscore")
    parser.add_argument("Gfa_file", help="Gfa file you want to extract a subgraph from")
    parser.add_argument("mode",choices=['contiguous', 'component'], help="selection mode : either find all 'contiguous' path starting from annotation file, or select all component including annotation")
    parser.add_argument("-o", help="output file name")
    parser.add_argument("-l",default="10000", help="max distance from annotation in bp (only for contiguous)")
    args = parser.parse_args()
    Annotation_file=args.Annotation_file
    List_nodes = list({line.rstrip().split("\t")[0] for index,line in enumerate(open(Annotation_file)) if index>0})
    List_nodes=[element for element in List_nodes if element]

    # List_nodes=list({"_".join(orf.split("\t")[0].split('_')[:-1]) for orf in open(Annotation_file)})
    # List_nodes=[element for element in List_nodes if element]

    Gfa_file=args.Gfa_file
    mode=args.mode
    if args.o :
        ouptut_file=args.o
    else :
        ouptut_file=Gfa_file.replace(".gfa","_"+Node_filename+".gfa")
    dist = int(args.l)
    graph_selection(List_nodes,Gfa_file,mode,dist,ouptut_file)
