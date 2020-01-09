# Metahood
Metahood is a pipeline entirely based on snakemake. 

**What does the  pipeline do :**
 - sample qualitycheck/trimming
- assemblies / co-assemblies
- binning (Concoct/Metabat2)
- de novo tree construction for mags
- diamond annotation and profiles
- output annotated orf graphs (derived from assembly graph), TO_FIX
- Strain resolution (Desman)

 **What we want to add :**
  - documentation
 - human Dna/contamination removal 
 - taxonomy profiling (CAT, kraken, ...) 
 - other options for  binning, e.g. maxbins2  
 - other bins assessment tools, e.g CheckM, Busco 
 - mags annotation and profiles
 
 **Overview of the rules workflows**
 This graph represent the binning part of the workflow strating from sample trimming.
![alt tag](./Binning.png)

###  How to install Metahood:
You can have a look at 
[conda_env.yaml](https://github.com/Sebastien-Raguideau/Metahood/blob/master/Conda_envs/conda_env.yaml)
, for an exhaustive list of all dependencies.  An already resolved environment is available in this same folder and can be used with this command line. 
```
cd path_to_repos/Metahood
conda env create -f Conda_envs/conda_env_MetaHood.yaml
```
You then need to activate the corresponding environment using : 

    conda activate MetaHood

**Fix CONCOCT install**
Unfortunately a bug still exist in the current conda package for concoct, you need to locate your conda installation and substitute it to `/var/lib/miniconda3/` in the following command :

    sed -i 's/original_data.values()/original_data.values/g' /var/lib/miniconda3/envs/MetaHood/bin/concoct_refine

**Databases** 
We rely on COG rpsblast database for MAG quality assessment. 

    wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz
    tar -xvf Cdd_LE.tar.gz

##  How to run Metahood:


    conda activate MetaHood
    path_to_repos/Metahood/start.py <output folder> --config <config file> -t <nb threads> -s <snakemake options> 


 ### Configuration file
 
The apparent lack of parameters is deceiving as all the complexity is hidden in a configuration file.  
[config.yaml](https://github.com/Sebastien-Raguideau/Metahood/blob/master/config.yaml)

 ------ Resssources ------
 **Threads** : Each task is allowed a maximum of 8 cores by default, you can change this value.
 **Percent_memory**: some tasks are memory intensive and some heuristic have been applied so that they don't use more than this number.
------ Path to data folder ------
**data**: Path to samples. Currently Metahood require for all samples to be stored in independant folders. The folder name will later define sample names in profiles. Samples must be at the format .fastq.gz and only 2 .fastq.gz can be present in any folder. 
------ Assembly parameters ------
**assembly**:
<p> 
**assembler** :  Currently only accept megahit.
**parameters**: any parameter you wish to pass to megahit
**per_sample**:  This line is optional, if you keep it, assembly and binning will also be carried per sample. You may specify a folder where to store these and also select the set of samples you want to have assembled, for instance : [per_sampleA|sampleA*] will create a per_sampleA directory and run a single sample approach on all samples starting with sampleA.
**groups**:  This line is optional, you may here specify any number of group of samples 
</p>
