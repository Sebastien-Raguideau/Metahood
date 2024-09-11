#!/usr/bin/env python3
import argparse
import numpy as np
from collections import defaultdict

def matrix_write_custom(matrix,file_name,col_names,row_names):
    ranks = ["root","Domain","Phylum","Class","Order","Genus","Family","Species"] 
    with open(file_name,"w") as handle:
        handle.write("%s\t%s\n"%("\t".join(ranks),"\t".join(col_names)))
        handle.writelines('%s\t%s\n'%("\t".join(row_names[index]),"\t".join(["{:.4g}".format(el) for el in line])) for index,line in enumerate(matrix))

def translate_sample(asmbl_path,sample):
    asmbl = asmbl_path.replace("/","_")
    if asmbl not in sample:
        sample = "%s_%s"%(asmbl,sample)
    return sample

def fill_taxa(taxa):
    return tuple(list(taxa)+(8-len(taxa))*["unknown"])



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("asm", help="list of assemblies")
    parser.add_argument("output", help="output")
    args = parser.parse_args()
    ASMBL = sorted((args.asm).split("-custom_sep-"))
    output = args.output

    get_covfile = lambda asmbl: "%s/profile/coverage_contigs.tsv"%(asmbl)
    get_bed = lambda asmbl: "%s/annotation/contigs.bed"%(asmbl)
    get_CAT = lambda asmbl: "%s/annotation/CAT_contigs_taxonomy.tsv"%(asmbl)

    # parse taxonomic annotation from CAT:
    spec_contigs = defaultdict(lambda:defaultdict(list))
    contig_spec = defaultdict(lambda:{})
    for asm in ASMBL:
        with open(get_CAT(asm)) as handle:
            _ = next(handle)
            for line in handle:
                sline = line.rstrip().split("\t")
                if len(sline)>3:
                    taxa = fill_taxa(sline[3].split(";"))
                else:
                    taxa = fill_taxa(["unknown"])
                contig_spec[asm][sline[0]] = taxa
                spec_contigs[asm][tuple(taxa)].append(sline[0])


    # get size of contigs in order to get total nuc per taxa
    contig_len = defaultdict(lambda:{})
    for asm in ASMBL:
        with open(get_bed(asm)) as handle:
            for line in handle:
                sline = line.rstrip().split("\t")
                contig_len[asm][sline[0]] = float(sline[2])

    # first create a matrix of the right dimension, mostly get merge diff entries in diff asmbl
    header_tot=[]
    for asmbl in ASMBL:
        with open(get_covfile(asmbl)) as handle:
            header = [translate_sample(asmbl,el) for el in next(handle).rstrip().split('\t')[1:]]
            header_tot+=header

    sorted_header = sorted(header_tot)
    taxa = {spec for asm,spec_to_cnt in spec_contigs.items() for spec in spec_to_cnt.keys()}
    sorted_annot = sorted(taxa)
    annot_to_index = {annot:index for index,annot in enumerate(sorted_annot)}
    annot_mat = np.zeros((len(sorted_annot),len(sorted_header)))

    # Then populate it by itering other going over assemblies
    for asmbl in ASMBL:
        with open(get_covfile(asmbl)) as handle:
            head_local = [translate_sample(asmbl,el) for el in next(handle).rstrip().split('\t')[1:]]
            col_ind = np.array([sorted_header.index(el) for el in head_local])
            for line in handle:
                splitline = line.rstrip().split("\t")
                row_index = annot_to_index[contig_spec[asmbl][splitline[0]]]
                annot_mat[row_index,col_ind] += np.array(list(map(float,splitline[1:])))*contig_len[asmbl][splitline[0]]

    # output matrix
    matrix_write_custom(annot_mat,output,sorted_header,sorted_annot)


