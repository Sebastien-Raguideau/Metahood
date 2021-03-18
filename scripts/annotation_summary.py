#!/usr/bin/env python3

from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
import argparse
import numpy as np
from collections importCounter
import glob

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", help="annotation")
    parser.add_argument("-h", help="best hits")
    parser.add_argument("-o", help="output")
    args = parser.parse_args()
    annotation=args.a
    BEST_HITS=args.h
    output=args.o

    # get annotation files for each mag
    folder=dirname(dirname(annotation[0]))
    annotations=BEST_HITS
    anno_to_files={anno:glob.glob("%s/*/*%s_best_hits.tsv"%(folder,anno)) for anno in annotations}
    # build annotation matrix for each annotation
    sorted_mags=sorted({file.split("/")[-2] for files in  anno_to_files.values() for file in files })
    folder_out=dirname(output[0])
    for anno, list_files in  anno_to_files.items() :
        annotation_count={basename(file).split('_%s'%anno)[0]:Counter([line.split('\t')[1] for line in open(file)]) for file in list_files }
        sorted_annotation=sorted({key for count in annotation_count.values() for key in count})
        annotation_to_index={key:index for index,key in enumerate(sorted_annotation)}
        matrix=np.zeros((len(sorted_mags),len(sorted_annotation)))
        for index_row,mag in enumerate(sorted_mags) :
            for key,count in annotation_count[mag].items() :
                index_col=annotation_to_index[key]
                matrix[index_row,index_col]= count
        with open("%s/%s_summary.tsv"%(folder_out,anno),"w") as handle : 
            handle.write("\t".join(["summary"]+sorted_annotation)+"\n")
            handle.write("\n".join(["\t".join([mag]+[str(int(nb)) for nb in matrix[index_row,:]]) for index_row,mag in enumerate(sorted_mags)])+"\n")

