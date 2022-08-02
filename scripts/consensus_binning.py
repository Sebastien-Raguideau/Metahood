#!/usr/bin/env python
# -*- coding: latin-1 -*-
import os
import argparse
import itertools
import numpy as np
from os.path import basename, dirname
from scipy.stats import pearsonr
from collections import Counter, defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp

# --------- basic algo -----------
# 1) focus on mags only, don't try to patch together incomplete bins or overly contaminated bins
# --> find contigs common between the 2 set of mags from the 2 assemblies
# 2) divide contigs in 2 situations:
#     2a) linked mags are the same thing/organism, just binning is different
#     2b) linked mags are different but got common contigs by chance
# ==> do they share at least 50% of their SCG. 50% is chosen because in the worst situation, to be MAG they need 75% completion with 0% contamination. In this case and unless they would still share 50% of their scgs. And we are assured that with that threshold we are not missing any pairing, while over 50% we could. Taking under 50% would also be meaningless since genuine non chimeric MAGs share at least 50%. Chimeric MAG: a bin falling in the MAG category from complementarity of multiple scg sets.
# 3) deal with contigs from case 2b), look at contigs coverage and mean cov of mags, assign contig to mag with least euclidian distance
# 4) Reassess quality, it is possible one of this bin lost it's scgs and droped from MAG to bin.
# 5) Deal with 2a) now that all other ambiguous contigs are sorted. compare pair of mags from both assembly and decide which is best from Completion -5*contamination criterion. 
# --> when this is done there is no more contig ambiguity even if some contigs may have been dropped.

# --------- more complexe situation -----------
# what happens if at 2a), it's more than 2 mags sharing at least 50% of their scg?
# what configuration may happens? 
# --> simple configuration, is 1 contaminated with 25% additional SCG sharing 50% of those with 2 non contaminated.
# --> worst case scenario everything is at a contamination of 25% and share 50% contigs with 2 mags from the other binner. Themselves with 25% contamination, sharing 50% of scg with 2 mags from the other binner.... 
# ----> solution just get a score per mag and remove iteratively the worst mag until no linkage remain.
# ----> there is a risk you will loose SCG and organism this way

def consensus(cluster_def_m2, cluster_def_c, profile_file, contig_to_len, contigs_to_scg,output):
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
    non_binary_relationship = set()
    for mag,linked_mag_scg in temp_smags.items() : 
        if linked_mag_scg=={} :
            continue
        if len(linked_mag_scg)>1:
            non_binary_relationship|={tuple(sorted([mag]+list(linked_mag_scg.keys())))}
        smags|={tuple(sorted((mag,mag_linked))) for mag_linked in linked_mag_scg.keys()}


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
    nb_dot = Counter([len(contig.split(".")) for contig in contigs])
    assert len(nb_dot)<=1,'nb of "." seems to be inconsistant in contigs names, this is an issue as there is no way to know how to recognise a split contig from concoct from a regular contigs.%s'%nb_dot.items()
    if len(nb_dot)==1:
        nb_dot=list(nb_dot.keys())[0]
        with open(profile_file) as handle :
            header=next(handle)
            contig_profile = {line.split(",")[0]:np.array(list(map(float,line.rstrip().split(",")[1:]))) for line in handle if ".".join(line.split(",")[0].split(".")[:nb_dot]) in contigs}
    else:
        with open(profile_file) as handle :
            header=next(handle)
            contig_profile = {line.split(",")[0]:np.array(list(map(float,line.rstrip().split(",")[1:]))) for line in handle if line.split(",")[0] in contigs}


    # we're reading concoct coverage file, it got splited contigs, we're going to take the means of the split's coverage
    contig_profile2 = defaultdict(list)
    for contig,profile in contig_profile.items() :
        if nb_dot!=0:
            contig = ".".join(contig.split(".")[:nb_dot])
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
    mags_to_delete = set(np.array(sorted_mags)[np.where((scg_tables==1).sum(1)<(0.75*36))])

    ## choose best representative from smags
    # define a criterion
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
        if not {mag1,mag2}&set(itertools.chain(*non_binary_relationship)):
            best_mag = get_best_mag(mag1, mag2, sorted_mags, mags_to_delete, contig_to_len, cluster_to_contigs)
            # delete the one which is not the best.... next time let's make it even more convoluted 
            mags_to_delete|={[mag1,mag2][mag1==best_mag]}

    # deal with non binary relationship: remove the worst of the serie, and then again....
    smags_remaining = lambda mtd:{el for el in smags if not mtd&set(el)}
    while smags_remaining(mags_to_delete):
        for case in non_binary_relationship:
            scores = [criterion(mag, sorted_mags, mags_to_delete) for mag in case]
            mags_to_delete |={case[min(range(len(scores)), key=scores.__getitem__)]}

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
        handle.writelines("%s,%s\n"%(contig,mag[0]) for contig,mag in cluster_def_final.items())

    # output the cluster definition in the same format as input : 
    with open("%s/mags_in_both_mapping.tsv"%dirname(output),'w') as handle : 
        handle.write("metabat2\tconcoct\n")
        handle.writelines("%s\t%s\n"%(*sorted(smag,key=lambda x:int(x[0]=="c")),) for smag in smags)


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
    consensus(cluster_def_m2, cluster_def_c, profile_file, contig_to_len, contigs_to_scg, output)





