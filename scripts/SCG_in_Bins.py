#!/usr/bin/env python
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import Counter, defaultdict
import argparse
import os


def get_contig_split(contig_bed, contig, scg_bed, orf) : 
    start,end = scg_bed[orf]
    if end<start : 
        (end,start) = (start,end)
    best_overlap = []
    for split_contig,start_split,end_split in contig_bed[contig] :
        if end_split<start_split : 
            (end_split,start_split) = (start_split,end_split)
        # test if included in split 
        is_start_in = (start_split<=start<=end_split)
        is_end_in = (start_split<=end<=end_split)
        # simplest case : 
        if is_end_in&is_start_in :
            return split_contig 
        # shit overlap case
        elif is_end_in|is_start_in :
            best_overlap.append([split_contig,min(end_split,end)-max(start_split,start)])
    return max(best_overlap,key=lambda x:x[1])[0]


def main(Bin_file, Fasta_file, C10K_bed, orf_bed, path, Table, LIST):
    # map contig + SCG name and sequence
    Dico_contig_SCGSeq = defaultdict(lambda: defaultdict(list))
    for header, seq in SimpleFastaParser(open(Fasta_file)):
        contig = "_".join(header.split(" ")[0].split('_')[:-1])
        SCG = header.split(" ")[1]
        Dico_contig_SCGSeq[contig][SCG].append([header, seq])
    
    # ---- deal with case where we have split contigs  ------
    # get split contigs bed infos 
    contig_bed = defaultdict(list)
    for line in open(C10K_bed) :
        contig,start,end,split_contig = line.rstrip().split("\t")
        if (len(split_contig.split("."))==1)|(contig not in Dico_contig_SCGSeq) :
            continue
        contig_bed[contig].append([split_contig,int(start),int(end)])
    # get SCG orf definitions 
    scg_bed = {line.rstrip().split("\t")[3]:list(map(int,line.rstrip().split("\t")[1:3])) for line in open(orf_bed) if line.rstrip().split("\t")[0] in Dico_contig_SCGSeq}
    # for each split contigs assign its SCG 
    Dico_contig_SCGSeq_splits = defaultdict(lambda: defaultdict(list))
    for contig, dict_scg in Dico_contig_SCGSeq.items() : 
        if contig not in contig_bed :
            continue
        for scg, list_header in dict_scg.items() :
            for header, seq in list_header :
                orf = header.split(" ")[0]
                split_contig = get_contig_split(contig_bed, contig, scg_bed, orf)
                Dico_contig_SCGSeq_splits[split_contig][scg].append([header, seq])
    Dico_contig_SCGSeq.update(Dico_contig_SCGSeq_splits)
    # map bins to SCG
    Dico_bins_SCG = defaultdict(lambda: defaultdict(list))
    Dico_bins_nbcontigs = defaultdict(int)
    with open(Bin_file) as handle : 
        _=next(handle)
        for line in handle:
            contig, Bin = line.rstrip().split(',')
            Dico_bins_nbcontigs[Bin] += 1
            if contig in Dico_contig_SCGSeq :
                for SCG, list_values in Dico_contig_SCGSeq[contig].items():
                    Dico_bins_SCG[Bin][SCG] += list_values
# --------------- SCG output for concoct refine------------------------------------------------------------
    if Table:
        List_SCG = sorted({key for dict in Dico_bins_SCG.values()
                           for key in dict.keys()})
        Dico_bin_Array = {bin_nb: [len(Dico_bins_SCG[bin_nb][SCG]) if Dico_bins_SCG[bin_nb]
                                   else 0 for SCG in List_SCG] for bin_nb in Dico_bins_nbcontigs}
        # add exception for consensus 
        def sort_bin_name(x) :
            if x[0].isalpha() :
                return int(x[1:])
            else :
                return int(x) 
        List_bins = sorted(Dico_bins_nbcontigs.keys(), key=sort_bin_name)

        SCG_table = [[Bin]+Dico_bin_Array[Bin] for Bin in List_bins]
        Handle = open(Table, "w")
        Handle.write(",".join(["Cluster"]+List_SCG)+"\n")
        Handle.write(
            "\n".join(map(lambda List: ','.join(map(str, List)), SCG_table)))
        Handle.close()

# -------------- create a folder by Mag with a folder by COG and their sequences inside -------------------
    if path:
        # which are 75% complete
        List_Mags = [Bin for Bin, List_contigs in Dico_bins_SCG.items() if sum(map(lambda x:x == 1, Counter([SCG for SCG, list_fasta in List_contigs.items() for header, seq in list_fasta]).values())) >= 0.75*36]
        # create a folder by Mag with a folder by COG and their sequences inside
        for Mag in List_Mags:
            List_contigs_SCG = []
            Mag_path = path+"Bin_"+Mag
            os.system("mkdir -p "+Mag_path)
            Handle = open(Mag_path+"/SCG.fna", "w")
            Handle.write("".join(map(lambda x: ">"+x[0]+"\n"+x[1]+"\n", [
                         Fasta for List_fasta in Dico_bins_SCG[Mag].values() for Fasta in List_fasta])))
            Handle.close()

# -------------- output a list of all Mags ----------------------------------------------------------------
    if LIST:
        # which are 75% complete
        List_Mags = [Bin for Bin, List_contigs in Dico_bins_SCG.items() if sum(map(lambda x:x == 1, Counter(
            [SCG for SCG, list_fasta in List_contigs.items() for header, seq in list_fasta]).values())) >= (0.75*36)]
        # create a folder by Mag with a folder by COG and their sequences inside
        Handle = open(LIST, "w")
        Handle.write("\n".join(List_Mags))
        Handle.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("Bin_file", help="Binning result file, csv, first column is the contig, second column is the bin")
    parser.add_argument("SCG_Fasta", help="fasta file of Orfs annotated as SCG")
    parser.add_argument("orf_bed", help="bed file of orfs definition : needed to handle concoct cut contigs")
    parser.add_argument("C10K_bed", help="bed file of cut contigs definition : needed to handle concoct cut contigs ")
    parser.add_argument("-f", help="output SCG sequences, one file by bin, takes the path to where you want to store the mags ")
    parser.add_argument("-l", help="output the list of Mags, takes names of the file as argument")
    parser.add_argument("-t", help="output SCG table, takes the name of a the table")
    args = parser.parse_args()
    Bin_file = args.Bin_file
    Fasta_file = args.SCG_Fasta
    C10K_bed = args.C10K_bed
    orf_bed = args.orf_bed
    path = ""
    Table = ""
    List = ""
    if args.f:
        path = args.f
    if args.t:
        Table = args.t
    if args.l:
        List = args.l
    main(Bin_file, Fasta_file, C10K_bed, orf_bed, path, Table, List)
