
# --------- Start of MagAnalysis -----------------------------------------------------------------

rule aggregate_SCG:
    input: "{group}/binning/folder_done"
    output: expand("{{group}}/binning/MagAnalysis/SCGs/{COG}.fna",COG=LIST_COGS) 
    params: "{group}/binning/"
    shell: """
     {SCRIPTS}/ExtractSCGs.py {params} {SCG_DATA}/scg_cogs_min0.97_max1.03_unique_genera.txt {params}MagAnalysis/
      """

rule extract_ref_SCGs:
    input:SCG_DATA+"/All_COGs.tar.gz"
    output: expand(SCG_DATA+"/All_{COG}.ffn",COG=LIST_COGS)
    shell:"tar -xf {input} -C {SCG_DATA}/ "

rule run_mafft:
    input: scg="{path}/SCGs/{COG}.fna",
           database_scg=SCG_DATA+"/All_{COG}.ffn"
    output: all_scg="{path}/AlignAll/{COG}_all.ffn",
            alignement="{path}/AlignAll/{COG}_all.gffn"
    log: "{path}/AlignAll/mafft_{COG}.log"
    threads:1000
    shell: """
    cat {input.scg} {input.database_scg} > {output.all_scg}
    mafft --thread {threads} {output.all_scg} > {output.alignement} 2>{log}
    """
rule trimal :
    input: "{path}/{COG}_all.gffn"
    output: "{path}/{COG}_al.gfa"
    shell: "trimal -in {input} -out {output} -gt 0.9 -cons 60"


rule generate_data_for_fastree:
    input: trimed=expand("{{path}}/{COG}_al.gfa",COG=LIST_COGS),
           alignement=expand("{{path}}/{COG}_all.gffn",COG=LIST_COGS)
    output: gfa="{path}/AlignAllR.gfa",
            names="{path}/Names.txt"
    shell:"""
    cat {input.alignement} | grep ">" | sed 's/_COG.*//' | sort | uniq | sed 's/>//g' > {output.names}
    {SCRIPTS}/MapTI.pl {SCG_DATA}/TaxaSpeciesR.txt < $({SCRIPTS}/CombineGenes.pl {output.names} {input.trimed}) > {output.gfa}
    """

rule Launch_FastTree:
    input: "{path}/AlignAllR.gfa"
    output: "{path}/AlignAllR.tree"
    log: "{path}/SelectR.out"
    shell:"""
    FastTreeMP -nt -gtr < {input} 2> {log} > {output}
    """



