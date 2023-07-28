#!/usr/bin/env python3
import argparse
import sys
import numpy as np
import pysam

from collections import defaultdict
from collections import Counter


def main(bam_file,out_stats,out_hists,maxlen):
    samfile = pysam.AlignmentFile(bam_file, "r")
    contigs_10K = {dicct["SN"]:dicct["LN"] for dicct in samfile.header.as_dict()["SQ"] if dicct["LN"]>=10000}
    with open(out_stats,"w") as handle_stat, open(out_hists,"w") as handle_hists:
        handle_hists.write("contig\tposisition\tcoverage\n")
        handle_hists.write("contig\tq50\t(q75-q25)/q50\n")
        for contig,length in contigs_10K.items():
            # hist
            pos_to_cnt = defaultdict(int)
            pos_to_cnt.update({col.reference_pos:col.nsegments for col in samfile.pileup(contig)})
            handle_hists.writelines("%s\t%s\t%s\n"%(contig,pos,pos_to_cnt[pos]) for pos in range(length))
            # stats
            depthArray = np.array(list(pos_to_cnt.values()))
            (lq,med,uq) = np.quantile(depthArray,[0.25,0.5,0.75])
            devF = (uq - lq)/med
            handle_stat.write("%s\t%s\t%s\n"%(contig,med,devF))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bam_file", help="bam file to analyse")
    parser.add_argument("cov_stat", help="output file, (q75-q25)/q50 for each contigs > 10000")
    parser.add_argument("cov_hists", help="output file, cov_histogram for each contigs > 10000")
    parser.add_argument("-l", help="contig size to use",default=10000)
    args = parser.parse_args()

    bam_file = args.bam_file
    out_stats = args.cov_stat
    out_hists = args.cov_hists
    maxlen = int(args.l)
    main(bam_file,out_stats,out_hists,maxlen)
