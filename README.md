# Metahood / snakemake tutorial

## Metahood
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
You can have a look at Conda_envs/conda_env.yaml, for an exhaustive list of all dependencies.  An already resolved environment is available in this same folder and can be used with this command line. 
```
conda env create -f Conda_envs/conda_env_MetaHood.yaml
```
You then need to activate the corresponding environment using : 

    conda activate MetaHood


###  How to run Metahood:


    conda activate MetaHood
    ~/repos/Metahood/start.py <output folder> --config <config file> -t <nb threads> -s <snakemake options> 


 **Configuration file**
 The apparent lack of parameters is hidding all the complexity in a configuration file in .yaml format.  
[https://github.com/Sebastien-Raguideau/Metahood/blob/master/config.yaml](https://github.com/Sebastien-Raguideau/Metahood/blob/master/config.yaml)

