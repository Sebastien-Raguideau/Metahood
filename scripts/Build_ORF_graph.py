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
    kmer_len=0
    if len(List_match)>1 :
        print("more than one Kmer length : "+"\t".join(List_match)) 
    else :
        try:
            kmer_len=int(List_match[0][:-1])
        except:
            kmer_len=0
    List_lines=[]
    List_edges=[]
    for Contig,Line in Dico_contig_lines.items() :
        # [S,name,sequence,LN:i:117,KC:i:85209]
        KC_contig=0
        Contig_len=len(Line.split("\t")[2])
        Match_KC=re.search("KC:i:(\d+)",Line)
        if Match_KC :
            KC_contig=int(Match_KC.groups()[0])
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
        (_,Contig1,sgn1,Contig2,sgn2,M,*args)=Line.rstrip().split("\t")
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
    Handle.write("\n".join(List_lines+List_edges))
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





