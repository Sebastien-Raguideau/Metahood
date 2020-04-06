#!/usr/bin/env python3
# -*- coding: latin-1 -*-
import os 
import argparse
import numpy as np
from os.path import basename
from collections import Counter,defaultdict
from subprocess import Popen, PIPE


def collate(cov_files,output) :
	sorted_files = sorted(cov_files)
	samples_names = sorted([".".join(basename(path).split('.')[:-2]) for path in sorted_files])
	sorted_feature = [line.rstrip().split("\t")[3] for line in open(sorted_files[0])]
	cov_matrix = np.zeros((len(sorted_feature),len(samples_names)))
	# the same bed file is used to generate all these coverage files. In consequence the features are ordered in the same way in each file
	# TOFIX : maybe not optimal way of storgin/writting this data. 
	for index_col,file in enumerate(sorted_files) :
		print(file)
		for index_row,line in enumerate(open(file)) :
			feature, cov = line.rstrip().split("\t")[3:]
			cov_matrix[index_row, index_col] = cov
	# output : 
	with open(output,"w") as handle : 
		handle.write("\t".join(["feature"]+samples_names)+"\n")
		for index,line in enumerate(cov_matrix) : 
			handle.write("\t".join([sorted_feature[index]]+list(map(str,line)))+"\n")


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-o", required=True, help="output file name")
	parser.add_argument("-l", nargs="+", required=True , help="list of coverage files")
	args = parser.parse_args()
	output = args.o
	cov_files = args.l
	collate(cov_files,output)





