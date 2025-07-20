from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp

rule semibin_split_bed:
    input: contig = "{group}/contigs/contigs.fa"
    output: temp("{group}/annotation/split_semibin2.bedtemp")
    run:
        # taken from split_contigs from semibin2
        with open(output[0],"w") as handle:
            for h, seq in sfp(open(input["contig"])):
                h = h.split(" ")[0]
                L = len(seq)
                if L < MIN_CONTIG_SIZE:
                    continue
                half = L // 2
                h1 = h + "_1"
                seq1 = seq[:half]
                h2 = h + "_2"
                seq2 = seq[half:]
                handle.write(f"{h}\t0\t{half}\t{h1}\n")
                handle.write(f"{h}\t{half}\t{L}\t{h2}\n")



rule semibin_cov_reformat:
    input: lambda w:["{group}/map/%s.split_semibin2.cov"%sample for sample in COBINNING_SAMPLES[w.group]]
    output: "{group}/binning/semibin2/output/cov.done"
    run:
        foldout = dirname(output[0])
        shell("touch {output}")
        for file in input:
            sample = basename(file).split(".split_semibin2")[0]
            out = f"{foldout}/semibin2_cov_{sample}.txt"
            with open(file) as handle_r, open(out,"w") as handle_w:
                for line in handle_r:
                    sline = line.rstrip().split("\t")
                    handle_w.write("%s\n"%"\t".join(sline[-2:]))

rule faster_semibin2:
    input: cov = "{group}/binning/semibin2/output/cov.done",
           contigs = "{group}/contigs/contigs.fa"
    output: "{group}/binning/semibin2/output/contig_bins.tsv"
    params: out = "{group}/binning/semibin2/output",
    threads: 32
    resources:
        slurm_partition = get_resource("partition",mult=4),
        mem_mb = get_resource("mem",mult=4)
    singularity: "docker://quay.io/biocontainers/semibin:2.1.0--pyhdfd78af_0"
    shell: "SemiBin2 single_easy_bin -i {input.contigs} -a {params.out}/semibin2_cov_* -o {params.out} -t {threads} --compression=none  && touch {output}"




# rule semibin2:
#     input: bams = semibin2_input,
#            contigs = "{group}/contigs/contigs.fa"
#     params: fold = "{group}/map",
#             out = "{group}/binning/semibin2"
#     output: "{group}/binning/semibin2/output/contig_bins.tsv"
#     threads: 32
#     resources:
#         slurm_partition = get_resource("partition"),
#         mem_mb = get_resource("mem")
#     singularity: "docker://quay.io/biocontainers/semibin:2.1.0--pyhdfd78af_0"
#     shell: "SemiBin2 single_easy_bin -i {input.contigs} -b {params.fold}/*_mapped_sorted.bam -o {params.out}/output -t {threads} --compression=none  && touch {output}"

rule semibin2_post_processing:
    input: "{path}/semibin2/output/contig_bins.tsv"
    output:"{path}/semibin2/clustering_semibin2.csv"
    run:
        with open(output[0],"w") as handle_w:
            with open(input[0]) as handle:
                header = next(handle)
                handle_w.write("contig_id,0\n")
                for line in handle:
                    contig,bin_nb = line.rstrip().split("\t")
                    if bin_nb == "-1":
                        continue
                    handle_w.write(f"{contig},{bin_nb}\n")


rule get_consensus_binning2 :
    # only support 2 binner as of now
    input : c_bin_def = "{path}/binning/consensus/concoct_vs_metabat2/clustering_consensus.csv",
            s_bin_def = "{path}/binning/semibin2/clustering_semibin2.csv",
            c_mag_list = "{path}/binning/consensus/concoct_vs_metabat2/consensus_MAG_list.txt",
            s_mag_list = "{path}/binning/semibin2/semibin2_MAG_list.txt",
            scg = "{path}/annotation/contigs_SCG.fna",
            contig_profiles = "{path}/binning/concoct/original_data_gt%s.csv"%MIN_CONTIG_SIZE,
            contig_bed = "{path}/annotation/contigs.bed"
    output : "{path}/binning/consensus/clustering_consensus.csv"
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb = get_resource("mem")
    shell :"""
    {SCRIPTS}/consensus_binning.py -c_bin_def {input.c_bin_def} -m_bin_def {input.s_bin_def} -c_mag_list {input.c_mag_list} -m_mag_list {input.s_mag_list} -scg {input.scg} -contig_profiles {input.contig_profiles} -contig_bed {input.contig_bed} -o {output}
    """
