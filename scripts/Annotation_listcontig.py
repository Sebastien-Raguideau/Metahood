#!/usr/bin/env python
# -*- coding: latin-1 -*-

import argparse
from  collections import defaultdict 

def main(File_m8,output,Just_list,Set_Annotation) :
	Dico_Annotation_ORF=defaultdict(list)
	header=next(open(File_m8)).split("\t")[0][0]!="k"
	if header :
		Handle=open(File_m8)
		Header=next(Handle)
	else :
		Handle=open(File_m8)
	for line in Handle :
		Splitline=line.rstrip().split("\t")
		ORF=Splitline[0]
		Annotation=Splitline[1]
		if Set_Annotation :
			if Annotation in Set_Annotation :
				Dico_Annotation_ORF[Annotation].append(ORF)
		else :
			Dico_Annotation_ORF[Annotation].append(ORF)
	Handle.close()
	List_Annotation=sorted(Dico_Annotation_ORF.keys())
	if Just_list :
		Handle=open(output,"w")
		Handle.write("\n".join(sorted([value for List in Dico_Annotation_ORF.values() for value in List])))
		Handle.close()
	else :	
		Handle=open(output,"w")
		Handle.write("\n".join(["\t".join([Annotation]+Dico_Annotation_ORF[Annotation]) for Annotation in List_Annotation]))
		Handle.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("Annotation_file", help="Should be a m8 file like tsv, where the first column correspond to the querry and second to the refference (annotation)") 
	parser.add_argument("Output", help="name of the output file" )
	parser.add_argument("-g", help="output just a list of genes to get each gene profile",action="store_true" )
	parser.add_argument("-r", help="Restrict to a list of annotation (KO/ARO/COG...) specified in a text file, one by line")
	args = parser.parse_args()
	File_m8=args.Annotation_file
	output=args.Output
	if args.g :
		Just_list=1
	else :
		Just_list=0
	Set_Annotation=set()
	if args.g :
		Set_Annotation={element.rstrip() for element in open(args.g)}
	main(File_m8,output,Just_list,Set_Annotation)
