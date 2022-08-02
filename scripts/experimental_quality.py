#!/usr/bin/env python3
import argparse
from collections import defaultdict
from os.path import basename

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('output', help='output file', required=True)
    parser.add_argument('markers', nargs='+', help='checkm marker outputs', required=True)
    args = parser.parse_args()
    MARKERS = args.marker
    OUTPUT = args.output

    header = next(open(MARKERS[0]))
    bin_to_best = defaultdict(lambda:[["0" for i in range(len(header.split("\t")))],"NA"])
    best = lambda bin1,bin2 : max([bin1,bin2],key=lambda x:float(x[0][5])-float(x[0][6]))
    for file in MARKERS:
        kind = basename(file).replace("checkm_","").split(".")[0]
        with open(file) as handle:
            header = next(handle)
            for line in handle:
                splitlines = line.rstrip().split("\t")
                bin_to_best[splitlines[0]] = best(bin_to_best[splitlines[0]],[splitlines,kind])
    sorted_bins = sorted([_bin for _bin in bin_to_best],key = lambda x:[bin_to_best[x][1],float(bin_to_best[x][0][5])])

    with open(OUTPUT,"w") as handle:
        handle.write("%s\t%s"%(type,header))
        handle.writelines("%s\t%s\n"%(bin_to_best[_bin][1],"\t".join(bin_to_best[_bin][0])) for _bin in sorted_bins)
