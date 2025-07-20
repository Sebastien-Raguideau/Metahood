#!/usr/bin/env python
# -*- coding: latin-1 -*-

import argparse
from  collections import defaultdict 

def main(m8_file, orf_set ,output, Set_Annotation):
	annotation_orfs = defaultdict(list)
	header = next(open(m8_file)).rstrip().split("\t")
	with open(m8_file) as handle:
		for line in handle:
			sline = line.rstrip().split("\t")
			if len(sline)<2:
				continue
			orf, annot = sline[:2]
			if orf not in orf_set:
				continue
			if (Set_Annotation!=set())&(annot not in Set_Annotation):
				continue
			annotation_orfs[annot].append(orf)
	List_Annotation = sorted(annotation_orfs.keys())
	with open(output,"w") as handle:
		handle.writelines("%s\t%s\n"%(annot,"\t".join(annotation_orfs[annot])) for annot in List_Annotation)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("Annotation_file", help="Should be a m8 file like tsv, where the first column correspond to the querry and second to the refference (annotation)") 
	parser.add_argument("orf_bed",help="bed file of orfs definition")
	parser.add_argument("Output", help="name of the output file" )
	parser.add_argument("-r", help="Restrict to a list of annotation (KO/ARO/COG...) specified in a text file, one by line")
	args = parser.parse_args()
	m8_file=args.Annotation_file
	output=args.Output
	Set_Annotation=set()
	if args.r :
		Set_Annotation={el.rstrip() for el in open(args.r)}
	# set of orfs 
	orf_set = {line.rstrip().split("\t")[3] for line in open(args.orf_bed)}
	main(m8_file, orf_set ,output, Set_Annotation)
