#!/usr/bin/env python3
import argparse
import numpy as np


def translate_sample(asmbl_path,sample):
    asmbl = asmbl_path.replace("/","_")
    if asmbl not in sample:
        sample = "%s_%s"%(asmbl,sample)
    return sample

def matrix_write(matrix,file_name,col_names,row_names):
    with open(file_name,"w") as handle:
        handle.write("/\t%s\n"%"\t".join(col_names))
        handle.writelines('%s\t%s\n'%(row_names[index],"\t".join(["{:.4g}".format(el) for el in line])) for index,line in enumerate(matrix))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("asm", help="list of assemblies")
    parser.add_argument("output", help="output")
    args = parser.parse_args()
    ASMBL = sorted((args.asm).split("-custom_sep-"))
    output = args.output

    header_tot=[]
    get_norm = lambda asmbl:"%s/profile/Normalisation.tsv"%(asmbl)
    for asmbl in ASMBL:
        with open(get_norm(asmbl)) as handle:
            header_tot+= [translate_sample(asmbl,el) for el in next(handle).rstrip().split('\t')[1:]]

    sorted_header = sorted(header_tot)

    # next step is to create a normalisation matrix
    norm = np.zeros((2,len(header_tot)))
    for asmbl in ASMBL:
        with open(get_norm(asmbl)) as handle:
            head_local = [translate_sample(asmbl,el) for el in next(handle).rstrip().split('\t')[1:]]
            col_ind = np.array([sorted_header.index(el) for el in head_local])
            for index,line in enumerate(handle):
                norm[index,col_ind] = list(map(float,line.rstrip().split('\t')[1:]))
    # output norm
    matrix_write(norm,output,sorted_header,["Nucleotides","median_scg"])