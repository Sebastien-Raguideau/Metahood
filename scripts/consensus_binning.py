#!/usr/bin/env python
# -*- coding: latin-1 -*-
import os
import argparse
import itertools
import numpy as np
from os.path import basename
from scipy.stats import pearsonr
from collections import Counter, defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp


def consensus(mags_m2, mags_c, cluster_def_m2, cluster_def_c, profile_file, contig_to_len, contigs_to_scg,output):
    # shared contigs amongs mags : 
    common_contigs =set(cluster_def_m2.keys())&set(cluster_def_c.keys())
    mag_shared = defaultdict(lambda :defaultdict(list))
    for contig in common_contigs :
        m2 = "m%s"%cluster_def_m2[contig]
        c = "c%s"%cluster_def_c[contig]
        mag_shared[m2][c].append(contig)
        mag_shared[c][m2].append(contig)

    # shared SCG among mags :
    contigs_to_scgnb = defaultdict(int)
    contigs_to_scgnb.update({contig:len(cogs) for contig,cogs in contigs_to_scg.items()})
    mag_shared_scg = {mag:{linked_mag:sum([contigs_to_scgnb[contig] for contig in contigs]) for linked_mag,contigs in linked_mag_contigs.items()} for mag,linked_mag_contigs in mag_shared.items()}

    # get mags found by both binner, let's call them  smags for shared mags : mag share more than 50%SCG 
    temp_smags = {mag:{linked_mag:nb_scg for linked_mag,nb_scg in linked_mag_scg.items() if nb_scg>0.5*36 }for mag,linked_mag_scg in mag_shared_scg.items()}
    smags = set()
    for mag,linked_mag_scg in temp_smags.items() : 
        if linked_mag_scg=={} :
            continue
        assert(len(linked_mag_scg)==1),"more than one potential unique association between mags from concoct and metabat2 : mag = %s, linked_mags = %s"%(mag,linked_mag) 
        smags|={tuple(sorted((mag,list(linked_mag_scg.keys())[0])))}

    # define list of ambiguous contigs shared between mags, for which an assignment is needed 
    # ignored pair in smags
    contig_to_mags = {}
    for mag1,linked_mag_contigs in mag_shared.items() :
        for mag2,contigs in linked_mag_contigs.items() :
            if tuple(sorted([mag1,mag2])) in smags :
                continue
            for contig in contigs : 
                contig_to_mags[contig]=[mag1,mag2]

    ### get median coordinate of all mags from concoct file :
    ## only go for contigs in ambiguous mags
    # get a unique mag definition :
    cluster_to_contigs = defaultdict(list)
    for contig,nb in cluster_def_m2.items() :
        cluster_to_contigs['m%s'%nb].append(contig)
    for contig,nb in cluster_def_c.items() :
        cluster_to_contigs['c%s'%nb].append(contig)
    ambiguous_mags = {mag for mags in contig_to_mags.values() for mag in mags}
    contigs = {contig for mag in ambiguous_mags for contig in cluster_to_contigs[mag]}
    with open(profile_file) as handle :
        header=next(handle)
        contig_profile = {line.split(",")[0]:np.array(list(map(float,line.rstrip().split(",")[1:]))) for line in handle if line.split(",")[0].split(".")[0] in contigs}

    # deal with splits, we're going to take the means of the splits
    contig_profile2 = defaultdict(list)
    for contig,profile in contig_profile.items() :
        contig = contig.split(".")[0]
        contig_profile2[contig].append(profile)
    contig_profile = {contig:np.array(profile).mean(0) for contig,profile in contig_profile2.items()}

    # get the median profile by mag, without ambiguous contigs
    mag_profile = {mag : np.median(np.array([contig_profile[contig] for contig in cluster_to_contigs[mag] if (contig not in contig_to_mags)&(contig in contig_profile)]),0) for mag in ambiguous_mags}

    # dispatch contigs to the most correlated mag, let's assign a contig to the best correlated mag
    for contig,[mag1,mag2] in  contig_to_mags.items() :
        worst_mag= [mag1,mag2][np.argmin([pearsonr(contig_profile[contig],mag_profile[mag1])[0],pearsonr(contig_profile[contig],mag_profile[mag2])[0]])]
        # delete contig from cluster definition of bad mag
        del cluster_to_contigs[worst_mag][cluster_to_contigs[worst_mag].index(contig)]

    ### reassess quality of mags 
    # build a scg table 
    cogs = sorted({ el for scgs in contigs_to_scg.values() for el in scgs})
    mag_to_scg = {mag:[contigs_to_scg[contig] for contig in mag_contigs] for mag,mag_contigs in cluster_to_contigs.items()}
    scg_tables = np.zeros((len(cluster_to_contigs),len(cogs)))
    sorted_mags = sorted(cluster_to_contigs.keys())
    for index, mag in enumerate(sorted_mags) :
        for contig in cluster_to_contigs[mag] :
            if contig in contigs_to_scg :
                for scg in contigs_to_scg[contig] :
                    scg_tables[index,cogs.index(scg)]+=1
    # ignore any mag  under the threshold : of at least 75% of mags in a unique copy
    mags_to_delete = set(np.array(sorted_mags)[np.where((scg_tables==1).sum(1)< 0.75*36)])

    ## choose best representative from smags
    # define a criterio
    def criterion(mag, sorted_mags, mags_to_delete) :
        if mag in mags_to_delete : 
            return -500
        contamination = sum(scg_tables[sorted_mags.index(mag),:]>1)/36.
        completion = 1-sum(scg_tables[sorted_mags.index(mag),:]==0)/36.
        score = completion - 5*contamination
        return score
    def get_best_mag(mag1, mag2, sorted_mags, mags_to_delete, contig_to_len, cluster_to_contigs) :
        score1 = criterion(mag1, sorted_mags, mags_to_delete)
        score2 = criterion(mag2, sorted_mags, mags_to_delete)
        if score1==score2 :
            len1 = sum([contig_to_len[contig] for contig in cluster_to_contigs[mag1]])
            len2 = sum([contig_to_len[contig] for contig in cluster_to_contigs[mag2]])
            # by default if they have the same length, mag1 will be chosen
            return [mag1,mag2][len1<len2]
        else :
            return [mag1,mag2][score1<score2]
    # get the list of mag to ignore :
    for mag1,mag2 in smags : 
        best_mag = get_best_mag(mag1, mag2, sorted_mags, mags_to_delete, contig_to_len, cluster_to_contigs)
        # delete the one which is not the best.... next time let's make it even more convoluted 
        mags_to_delete|={[mag1,mag2][mag1==best_mag]}

    # create the final cluster definition 
    cluster_to_contigs_final = {mag:contigs for mag,contigs in cluster_to_contigs.items() if mag not in mags_to_delete}

    # check the process did not fail : 
    cluster_def_final = defaultdict(list)
    for mag,contigs in cluster_to_contigs_final.items() :
        for contig in contigs : 
            cluster_def_final[contig].append(mag)
    assert(sum([len(el)==1 for el in cluster_def_final.values()])==len(cluster_def_final)),"some contigs are still ambiguous"

    # output the cluster definition in the same format as input : 
    with open(output,'w') as handle : 
        handle.write("contig_id,0\n")
        handle.write("\n".join(["%s,%s"%(contig,mag[0]) for contig,mag in cluster_def_final.items()])+"\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c_bin_def", required=True, help="concoct bin definition : csv, first column = contig, second = cluster")
    parser.add_argument("-m_bin_def", required=True, help="metabat2 bin definition : csv, first column = contig, second = cluster")
    parser.add_argument("-c_mag_list", required=True, help="concoct mag number, one number by line")
    parser.add_argument("-m_mag_list", required=True, help="metabat2 mag number, one number by line")
    parser.add_argument("-scg", required=True, help="SCG annotation of contigs in the for of a .fna. Each orf header must be space delimited with element : orf,cog,strand")
    parser.add_argument("-contig_profiles", required=True, help="concoct generated contig information (tetranucleotide profile and coverage), ex : original_data_gt1000.csv")
    parser.add_argument("-contig_bed", required=True, help="bogus bed used for computing contigs coverages in bedtool, the feature is the contig itself")
    parser.add_argument("-o", required=True, help="output, unambiguous bin definition")


    # get args 
    args = parser.parse_args()
    # get mags definition : 
    mags_m2 = {line.rstrip() for line in open(args.m_mag_list)}
    mags_c = {line.rstrip() for line in open(args.c_mag_list)}
    # get cluster definition : 
    cluster_def_m2 = {line.rstrip().split(",")[0]:line.rstrip().split(",")[1] for line in open(args.m_bin_def) if line.rstrip().split(",")[1] in mags_m2}
    cluster_def_c = {line.rstrip().split(",")[0]:line.rstrip().split(",")[1] for line in open(args.c_bin_def) if (line.rstrip().split(",")[1] in mags_c)&(line.rstrip().split(",")[0]!="contig_id")}
    # contig info :
    profile_file = args.contig_profiles 
    # get contig size : 
    contig_to_len = {line.rstrip().split("\t")[0]:int(line.rstrip().split("\t")[2]) for line in open(args.contig_bed)}
    # get contig SCG
    contigs_to_scg = defaultdict(list)
    for header,seq in sfp(open(args.scg)) :
        orf,cog,strand = header.rstrip().split(" ")
        contig = "_".join(orf.split('_')[:-1])
        contigs_to_scg[contig].append(cog)
    # output 
    output = args.o

    #main 
    consensus(mags_m2, mags_c, cluster_def_m2, cluster_def_c, profile_file, contig_to_len, contigs_to_scg, output)





