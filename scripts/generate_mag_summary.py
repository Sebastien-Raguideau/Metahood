#!/usr/bin/env python3
import glob
import argparse
import numpy as np
from os.path import dirname,basename
from collections import defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp

# sort mag per RED from smaller to bigger 
def red_sorting(x):
    red = mag_to_red[mag_to_dmag[x]]
    if isinstance(red, float):
        return red
    else:
        return 1

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
    parser.add_argument("ar", help="list of assemblies")
    parser.add_argument("bac", help="output")
    parser.add_argument("drep", help="output")
    parser.add_argument("scheme", help="output")
    parser.add_argument("asm", help="output")
    parser.add_argument("output", help="output")    

    args = parser.parse_args()
    path_drep = dirname(args.drep)
    bac_gtdb = args.bac
    ar_gtdb = args.ar
    ASMBL = sorted((args.asm).split("-custom_sep-"))
    mag_naming_scheme = args.scheme
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
    PATH = dirname(dirname(output))
    mag_to_length = {basename(file).replace(".fa",""):sum([len(seq) for _,seq in sfp(open(file))]) for file in glob.glob("%s/all_mags/*/*.fa"%PATH)}

    # get mag name
    mag_path = {line.rstrip().split("\t")[0]:line.rstrip().split("\t")[2] for index,line in enumerate(open(mag_naming_scheme)) if index>0}

    # get quality/contamination
    mag_to_qual = {}
    for asmbl_path in ASMBL:
        with open("%s/binning/consensus/consensus_SCG_table.csv"%asmbl_path) as handle:
            _ = next(handle)
            mag_old_to_new = {basename(path):mag for mag,path in mag_path.items() if asmbl_path in path}
            for line in handle:
                splitline = line.rstrip().split(",")
                mag_name = mag_old_to_new["Bin_%s.fa"%splitline[0]]
                mag_to_qual[mag_name] = list(map(int,splitline[1:]))

    mag_to_completion = {mag:sum(np.array(mag_to_qual[mag])>0)/float(len(mag_to_qual[mag])) for mag in mag_to_dmag}
    mag_to_contamination = {mag:sum(np.array(mag_to_qual[mag])>1)/float(len(mag_to_qual[mag])) for mag in mag_to_dmag}

    # get mag/dmag mapping
    mag_to_dmag,dmag_to_mags = get_mag_dmag("%s/drep/data_tables"%PATH)

    # dmag prevalence
    dmag_to_prevalence = {dmag.replace(".fa",""):len({mag.split("_Bin")[0] for mag in mags}) for dmag,mags in dmag_to_mags.items()}

    # mag to assembly
    get_asm = lambda mag:"/".join(mag_path[mag].split("/")[:-4])

    # sort by red
    sorted_mags = sorted([mag for mag in mag_to_completion], key=red_sorting)

    # output all infos :
    with open(output,'w') as handle:
        handle.write("mag\tassembly\tdmag\tprevalence\tMAG_length\tcompletion\tcontamination\tRED\tDomain\tPhylum\tClass\tOrder\tfamilly\tgenus\tspecies\n")
        handle.writelines("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(mag,get_asm(mag),mag_to_dmag[mag],dmag_to_prevalence[mag_to_dmag[mag]],mag_to_length[mag],mag_to_completion[mag],mag_to_contamination[mag],mag_to_red[mag_to_dmag[mag]],"\t".join(mag_to_taxa[mag_to_dmag[mag]])) for mag in sorted_mags)
