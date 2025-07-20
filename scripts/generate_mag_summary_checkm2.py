#!/usr/bin/env python3
import os
import glob
import argparse
import numpy as np
from functools import partial
from os.path import dirname,basename
from collections import defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp


# sort mag per RED from smaller to bigger 
def red_sorting(x,mag_to_red,mag_to_dmag):
    red = mag_to_red[mag_to_dmag[x]]
    if isinstance(red, float):
        return red
    else:
        return 1


# mag_length = sorted(mag_to_qual.items(),key=lambda x:-int(x[1][1]))
# mag_length2 = [(mag,q[0],q[1],test(q)) for mag,q in mag_length]

# def test(q):
#     [length,N50,gc,cd,comp,cont] = list(map(float,q[:-1]))
#     return 100*N50/length

def score(qual):
    [length,N50,gc,cd,comp,cont] = list(map(float,qual[:-1]))
    # 5% contamination ~ same impact as circular
    score = comp - 5*cont + 25*N50/length
    return score

def get_mag_dmag(path):
    mag_to_cluster = {line.rstrip().split(",")[0]:line.rstrip().split(",")[1] for index,line in enumerate(open("%s/Cdb.csv"%path)) if index>0}
    header = next(open("%s/Wdb.csv"%path)).split(",")
    cl_ind = header.index("cluster")
    g_ind = header.index("genome")
    cluster_to_dmag = {line.rstrip().split(",")[cl_ind]:line.rstrip().split(",")[g_ind] for line in open("%s/Wdb.csv"%path)}
    mag_to_dmag = {mag.replace(".fa",""):cluster_to_dmag[cluster].replace(".fa","") for mag,cluster in mag_to_cluster.items()}
    dmag_to_mags = defaultdict(list)
    for mag,dmag in mag_to_dmag.items():
        dmag_to_mags[dmag].append(mag)
    return mag_to_dmag,dmag_to_mags

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ar", help="")
    parser.add_argument("bac", help="")
    parser.add_argument("drep", help="")
    parser.add_argument("qual", help="")    
    parser.add_argument("output", help="")


    args = parser.parse_args()
    path_drep = dirname(args.drep)
    bac_gtdb = args.bac
    ar_gtdb = args.ar
    qual_file = args.qual
    output = args.output

    # get mag/dmag relationship
    mag_to_dmag,dmag_to_mags = get_mag_dmag(path_drep)

    # extract gtdb taxonomy
    mag_to_taxa = {line.rstrip().split("\t")[0]:[taxa.split("__")[1] for taxa in line.rstrip().split("\t")[1].split(';')] for index,line in enumerate(open(bac_gtdb)) if index>0}
    mag_to_taxa.update({line.rstrip().split("\t")[0]:[taxa.split("__")[1] for taxa in line.rstrip().split("\t")[1].split(';')] for index,line in enumerate(open(ar_gtdb)) if index>0})
    mag_to_taxa = {mag:[tax if tax else "/" for tax in taxa] for mag,taxa in mag_to_taxa.items()}
    mag_to_taxa2 = defaultdict(lambda:7*["/"])
    mag_to_taxa2.update(mag_to_taxa)
    mag_to_taxa = mag_to_taxa2

    # get mag red
    mag_to_red = defaultdict(lambda:"/")
    mag_to_red.update({line.rstrip().split("\t")[0]:float(line.rstrip().split("\t")[-2]) if line.rstrip().split('\t')[-2]!="N/A" else 1 for index,line in enumerate(open(bac_gtdb)) if index>0})
    mag_to_red.update({line.rstrip().split("\t")[0]:float(line.rstrip().split("\t")[-2]) if line.rstrip().split('\t')[-2]!="N/A" else 1 for index,line in enumerate(open(ar_gtdb)) if index>0})

    # get mag length
    mag_to_qual = {}
    with open(qual_file) as handle:
        header = next(handle)
        for line in handle:
            mag, comp, cont, _, _, cd, N50, _, length, gc, _, _,_,_ = line.rstrip().split("\t")
            cp = float(comp)
            ct = float(cont)
            quality=""
            if (ct>10):
                quality = "contaminated(>10)"                
            if (cp<50)|(ct>25):
                quality = "trash(<50/>25)"
            if (cp>=50)&(ct<=10):
                quality = "Low-quality(50/10)"
            if (cp>=75)&(ct<=10):
                quality = "Medium-quality(75/10)"
            if (cp>=90)&(ct<=5):
                quality = "High-quality(90/5)"
            if quality=="":
                to_fix.append([cp,ct])
            mag_to_qual[mag] = [length,N50,gc,cd,comp,cont,quality]


    # score mags and rechoose dmag
    updated_dmag_to_mags = {}
    updated_mag_to_dmag = {}
    for dmag,mags in dmag_to_mags.items():
        if dmag in mag_to_qual:
            scored_mags = sorted(mags,key=lambda x:-score(mag_to_qual[x]))
            new_dmag = scored_mags[0]
            updated_mag_to_dmag.update({m:new_dmag for m in mags})
            updated_dmag_to_mags[new_dmag] = mags
            mag_to_taxa[new_dmag] = mag_to_taxa[dmag]

    # be sure not to use old values
    del mag_to_dmag
    del dmag_to_mags

    # dmag prevalence
    dmag_to_prevalence = {dmag.replace(".fa",""):len({mag.split("_Bin")[0] for mag in mags}) for dmag,mags in updated_dmag_to_mags.items()}

    # mag to assembly
    get_asm = lambda mag:mag.split("_Bin")[0]

    # sort by red
    sorted_mags = sorted([mag for mag in mag_to_qual], key=partial(red_sorting,mag_to_red=mag_to_red,mag_to_dmag=updated_mag_to_dmag))

    # output all infos :
    with open(output,'w') as handle:
        handle.write("mag\tassembly\tdmag\tprevalence\tMAG_length\tN50\tGC\tcoding_density\tcompletion\tcontamination\tquality\tRED\tDomain\tPhylum\tClass\tOrder\tfamilly\tgenus\tspecies\n")
        handle.writelines("%s\n"%"\t".join([mag,get_asm(mag),updated_mag_to_dmag[mag],str(dmag_to_prevalence[updated_mag_to_dmag[mag]])]+mag_to_qual[mag]+[str(mag_to_red[updated_mag_to_dmag[mag]])]+ mag_to_taxa[updated_mag_to_dmag[mag]]) for mag in sorted_mags)

