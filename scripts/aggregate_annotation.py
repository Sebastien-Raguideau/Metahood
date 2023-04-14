#!/usr/bin/env python3
import argparse
import numpy as np

def matrix_write(matrix,file_name,col_names,row_names):
    with open(file_name,"w") as handle:
        handle.write("/\t%s\n"%"\t".join(col_names))
        handle.writelines('%s\t%s\n'%(row_names[index],"\t".join(["{:.4g}".format(el) for el in line])) for index,line in enumerate(matrix))

def translate_sample(asmbl_path,sample):
    asmbl = asmbl_path.replace("/","_")
    if asmbl not in sample:
        sample = "%s_%s"%(asmbl,sample)
    return sample


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("annot", help="type of annotation")
    parser.add_argument("asm", help="list of assemblies")
    parser.add_argument("output", help="output")
    args = parser.parse_args()
    annot = args.annot
    ASMBL = sorted((args.asm).split("-custom_sep-"))
    output = args.output


    get_covfile = lambda asmbl: "%s/profile/cov_%s.tsv"%(asmbl,annot)

    # first create a matrix of the right dimension, mostly get merge diff entries in diff asmbl
    header_tot=[]
    annotation = set()
    for asmbl in ASMBL:
        with open(get_covfile(asmbl)) as handle:
            header = [translate_sample(asmbl,el) for el in next(handle).rstrip().split('\t')[1:]]
            header_tot+=header
            annotation|={line.rstrip().split("\t")[0] for line in handle}

    sorted_header = sorted(header_tot)
    sorted_annot = sorted(annotation)
    annot_to_index = {annot:index for index,annot in enumerate(sorted_annot)}
    annot_mat = np.zeros((len(sorted_annot),len(sorted_header)))

    # Then populate it by itering other going over assemblies
    for asmbl in ASMBL:
        with open(get_covfile(asmbl)) as handle:
            head_local = [translate_sample(asmbl,el) for el in next(handle).rstrip().split('\t')[1:]]
            col_ind = np.array([sorted_header.index(el) for el in head_local])
            for line in handle:
                splitline = line.rstrip().split("\t")
                row_index = annot_to_index[splitline[0]]
                annot_mat[row_index,col_ind] +=list(map(float,splitline[1:]))

    # output matrix
    matrix_write(annot_mat,output,sorted_header,sorted_annot)
