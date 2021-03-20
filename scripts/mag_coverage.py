#!/usr/bin/env python3

import argparse
import numpy as np
from collections import defaultdict

def get_cluster_def(file,mags):
    cluster_def=defaultdict(list)
    contig_to_cluster={}
    for line in open(file):
        contig,cluster=line.rstrip().split(",")
        if cluster in mags:
            cluster_def[cluster].append(contig)
            contig_to_cluster[contig]=cluster
    return cluster_def,contig_to_cluster

def get_contig_cov(file,contigs):
    contig_cov={}
    with open(file) as handle :
        header=next(handle)
        for line in handle:
            splitline=line.rstrip().split("\t")
            contig=splitline[0]
            if contig in contigs:
                contig_cov[contig]=np.array(list(map(float,splitline[1:])))
    return header,contig_cov


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", help="mags list")
    parser.add_argument("-c", help="cluster csv")
    parser.add_argument("-t", help="cov tsv")
    parser.add_argument("-l", help="len contigs bed")
    parser.add_argument("-n", help="nb nuc tsv")
    parser.add_argument("-v", help="mag coverage")
    parser.add_argument("-p", help="mag map")
    args = parser.parse_args()

    mags_list=args.m
    cluster=args.c
    cov=args.t
    length=args.l
    nb_nuc=args.n
    mag_cov=args.v
    mag_map=args.p


    # get the sets of bins, considered mags
    mags={line.rstrip() for line in open(mags_list)}
    # get cluster definition
    cluster_def,contig_to_cluster=get_cluster_def(cluster,mags)
    # get contigs coverage
    header,contig_cov=get_contig_cov(cov,contig_to_cluster)
    # get contigs lengths
    contig_len={line.rstrip().split('\t')[0]:float(line.rstrip().split('\t')[2]) for line in open(length) if line.rstrip().split('\t')[0] in contig_to_cluster}
    # get nb reads :
    with open(nb_nuc) as handle:
        header_norm = next(handle).rstrip().split("\t")[1:]
        samples = header.rstrip().split('\t')[1:]
        nuc_nb=list(map(float,next(handle).rstrip().split("\t")[1:]))
        if samples != header_norm :
            reorder = [header_norm.index(elmt) for elmt in samples]
            header_norm = np.array(header_norm)[reorder]
            nuc_nb = np.array(nuc_nb)[reorder]
        assert samples == header_norm,"error, header are differents : \n%s\n%s"%(samples, header_norm)
     
    # get mag nucleotides
    handle_cov=open(mag_cov,"w")
    handle_map=open(mag_map,"w")
    handle_map.write("mag\t"+"\t".join(samples)+"\n")
    handle_cov.write("mag\t"+"\t".join(samples)+"\n")
    count = 0
    for cluster,list_contigs in cluster_def.items():
        if cluster =="0":
            print("Bin_%s\t"%cluster,count)
        tot_len=0
        tot_nuc=np.zeros(len(samples))
        for contig in list_contigs:
            if contig=="contig_id":
                continue
            tot_len+=contig_len[contig]
            tot_nuc+=contig_cov[contig]*contig_len[contig]
        handle_cov.write("Bin_%s\t"%cluster+"\t".join(list(map(str,tot_nuc/tot_len)))+"\n")
        handle_map.write("Bin_%s\t"%cluster+"\t".join(list(map(str,tot_nuc/nuc_nb)))+"\n")
        count+=1
    handle_cov.close()
    handle_map.close()

