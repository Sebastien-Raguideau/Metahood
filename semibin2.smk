
rule semibin2:
    input:  bams = lambda w : ["{path}/{group}/map/%s_mapped_sorted.bam"%basename(sample) for sample in  GROUPS[w.group]],
           contigs = "{path}/{group}/contigs/contigs.fa",
    params: fold = "{path}/{group}/binning/semibin2"
    output: "{path}/{group}/binning/semibin2/output/contig_bins.tsv"
    threads: 32
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb = get_resource("mem")
    singularity: "docker://quay.io/biocontainers/semibin:2.1.0--pyhdfd78af_0"
    shell: "SemiBin2 single_easy_bin -i {input.contigs} -b {params.fold}/bam/*.bam -o {params.fold}/output --sequencing-type=long_read -t {threads} --compression=none  && touch {output}"

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
