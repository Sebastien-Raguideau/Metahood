#!/usr/bin/env python3
import argparse
from collections import defaultdict,Counter


def get_contigs_ORF_list(Gfa_file) :
	List_ORFs=[]
	Handle=open(Gfa_file)
	for Line in Handle : 
		Line=Line.rstrip().split("\t")
		if Line[0]=="S" : 
			# [S,name,sequence,LN:i:117,KC:i:85209,CL:z:#5cb9be]
			name=Line[1]
			List_ORFs.append(name)
		if Line[0]=="L" :
			break
	Handle.close()
	return List_ORFs


def main(Gfa_file,Dico_annotation_color,output) :
	List_ordered_orfs=get_contigs_ORF_list(Gfa_file)
	dico_orf_annotation=defaultdict(list)
	Header_global=next(open(list(Dico_annotation_color.keys())[0])).rstrip().split('\t')[:7]
	Dico_index_header={}
	Set_orfs_annotation=set([])
	for index,file in enumerate(Dico_annotation_color) : 
		Handle=open(file)
		Header=next(Handle).rstrip().split('\t')
		Dico_index_header[index]=Header[7:]
		Header_global+=Header[7:]
	Dico_index_Padding={index:["" for nb in range(len(header))] for index,header in Dico_index_header.items()}
	List_index_dict_color=[]
	for index,(file,color) in enumerate(Dico_annotation_color.items()) :
		List_index_dict_color.append([index,{line.rstrip().split('\t')[0]:line.rstrip().split('\t') for line in open(file)},color])
	List_results=[Header_global+['Colour']]
	for ORF in List_ordered_orfs :
		color="#696969"
		Line_result=[[ORF]+["" for nb in range(6)]]+[[] for nb in List_index_dict_color]+[[color]]
		for index,Dict,color in List_index_dict_color :
			if ORF in Dict :
				Base_info=Dict[ORF][:7]
				Additional_info=Dict[ORF][7:]
				Line_result[0]=Base_info
				Line_result[1+index]=map(lambda x:x.replace(',','|'),Additional_info)
				Line_result[-1]=[color]
			else :
				Line_result[1+index]=Dico_index_Padding[index]
		True_line_result=[]
		for element in Line_result :
			True_line_result+=element
		List_results.append(True_line_result)
	Handle=open(output,"w")
	Handle.write("\n".join([",".join(Line) for Line in List_results]))
	Handle.close()


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("Gfa_file", help="gfa file for which we're adding color/annotation")
	parser.add_argument("colors", help="tsv file, first columns are annotation files, second columns are associated colors")
	parser.add_argument("output", help="output path/name")
	args = parser.parse_args()
	Gfa_file=args.Gfa_file
	Dico_annotation_color={line.rstrip().split("\t")[0]:line.rstrip().split("\t")[1] for line in open(args.colors)}
	output=args.output
	main(Gfa_file,Dico_annotation_color,output)

