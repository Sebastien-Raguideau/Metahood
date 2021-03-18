#!/usr/bin/env python3

from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", help="GROUPS")
    parser.add_argument("-o", help="output folder")
    args = parser.parse_args()
    GROUPS=args.g
    output_folder=args.o

    map_mag_to_bin = [["mag_name","assembly_nb","mag_nb","assembly_name","bin_name"]]
    mag_nb = 0
    for index_group, group in enumerate(GROUPS) : 
        path = "%s/binning/consensus"%group
        file = "%s/consensus_MAG_list.txt"%path
        bins_path = "%s/bins"%path
        for line in open(file) : 
            bin_name = line.rstrip()
            mag_nb+=1
            name = "a%s_m%s.fa"%(index_group,mag_nb)
            map_mag_to_bin.append([name.split(".fa")[0],index_group,mag_nb,group,"Bin_%s"%bin_name])
            with open("%s/%s"%(output_folder,name),"w") as handle :
                for header,seq in sfp(open("%s/Bin_%s.fa"%(bins_path,bin_name))) :
                    new_header = "assembly_%s_%s"%(index_group,header)
                    handle.write(">%s\n%s\n"%(new_header,seq))
    with open(output["map"],"w") as handle :
        handle.write("\n".join(["\t".join(list(map(str,line))) for line in map_mag_to_bin])+"\n")
