#!/bin/bash

container_list=("pythonenv" "bandage" "bedtools" "blast" "bwasamtools" "cat" "concoct" "desman" "diamond" "drep" "fastqc" "fasttree" "gtdbtk" "kofamscan" "krakenuniq" "mafft" "megahit" "metabat2" "multiqc" "prodigal" "trimal" "trim_galore")

mkdir -p builds

for item in ${container_list[@]}; do
    sudo singularity build builds/${item}.sif container_recipes/singularity/Singularity.${item}
done
