from __future__ import print_function
from __future__ import division

try:
    from future_builtins import zip
except:
    pass

from os.path import basename, dirname, realpath, isfile
from collections import defaultdict, Counter
from subprocess import Popen, PIPE
from functools import partial
import os.path
import glob
import re
import os


default_values = {
    "binning":{"concoct":{"contig_size" : 1000,"execution" : 1,"max_bin_nb" : 2000},"metabat2":{"execution" : 1,"contig_size":1500},"ssa_unique_sample":False,"cobinning_samples":["*"]},
    "mag":["native"],
    "threads":8,
    "assembly":    {"assembler": "megahit","groups": {},"parameters":"" },
    "annotation": {'diamond':{},"cat_db":"","cat_path":"","kraken_db":"","kofamscan":{"profiles":"","ko_list":""},"virsorter":"","plasmidnet_install":"","genomad_db":""},
    "graph":{"List_graphs":{}},
    "filtering":"",
    "Percent_memory":0.5,
    "task_memory":200,
    "maganalysis":0,
    "desman":{"execution":0, "nb_haplotypes": 10,"nb_repeat": 5,"min_cov": 1,"scripts":""},
    "slurm_partitions":{"":{"name":"","min_mem":"","max_mem":"","min_threads":"","max_threads":""}},
    "nb_map":20
}



# ---- neat regex matching of files --------
def extended_glob(pattern):
    process = Popen(['bash -O extglob -c " ls -d '+pattern+'/ "'], stdout=PIPE, stderr=PIPE,shell=True)
    List_path=[element[:-1] for element in process.communicate()[0].decode("utf-8").split("\n") if element]
    return [path for path in List_path if basename(path)!="multiqc_data"]


# Taken from http://stackoverflow.com/questions/36831998/how-to-fill-default-parameters-in-yaml-file-using-python
def setdefault_recursively(tgt, default = default_values):
    for k in default:
        if isinstance(default[k], dict): # if the current item is a dict,
            # expand it recursively
            setdefault_recursively(tgt.setdefault(k, {}), default[k])
        else:
            # ... otherwise simply set a default value if it's not set before
            tgt.setdefault(k, default[k])

def fill_default_values(config):
    local_dir = config.get("LOCAL_DIR")
    if local_dir:
        default_values["scripts"] = os.path.join(local_dir, "scripts")
        default_values["scg_data"] = os.path.join(local_dir, "scg_data")
        default_values["conda_env"] = os.path.join(local_dir, "conda_envs")
        default_values["gtdb"] = os.path.join(local_dir, "gtdb")
    setdefault_recursively(config)

def sample_name(fullname):
    return os.path.splitext(os.path.basename(fullname))[0]

# FASTA_EXTS = {".fastq", ".fastq.gz"}  # only 2 extension are valid. 
FASTA_EXTS = {".fastq.gz",".fastq",".fq.gz",".fq",".fa",".fasta",".fasta.gz",".fa.gz"}  # only extension valid. 

def gather_paths(path):
    for filename in os.listdir(path):
        name = os.path.basename(filename)
        for ext in FASTA_EXTS:
            if not name.endswith(ext):
                continue
            if "trimmed" in name : 
                continue
            if "Filtered" in name :
                continue
            yield os.path.join(path, filename)

def detect_reads(dir):
    return sorted(list(gather_paths(dir)))

def get_extension(file) :
    for ext in FASTA_EXTS:
        if file.endswith(ext):
            return ext 

def replace_extensions(sample,FILTER):
    if FILTER:
        sample = "%s/Filtered_%s"%(dirname(sample),basename(sample))
    ext = get_extension(sample)
    if ext in {".fastq.gz",".fastq",".fq.gz",".fq"} :
        return sample.replace(ext,"_trimmed%s"%ext)
    else :
        return sample

def samples_yaml():
    libs = []
    for s in SAMPLES:
        info = {}
        info["left reads"] = [SAMPLE_READS[s][0]]
        info["right reads"] = [SAMPLE_READS[s][1]]
        info["type"] = "paired-end"
        info["orientation"] = "fr"
        libs.append(info)
    return yaml.dump(libs, default_style='"', default_flow_style=False)


#copied from http://stackoverflow.com/questions/431684/how-do-i-cd-in-python/13197763#13197763
class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


        
# from doc: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html
# from https://stackoverflow.com/questions/50891407/snakemake-how-to-dynamically-set-memory-resource-based-on-input-file-size

def get_resource_real(wildcards, input, threads, attempt, SLURM_PARTITIONS="", mode="", mult=2, min_size=10000):
    # return either partition, mem or threads, theses need to be done together
    # because if mem/cpu/threads is over partition definition, then mem/threads needs to change
    # since there can be setting for minimum mem/threads in some partition definition
    # also max threads needs 
    def return_result(mem,partition,threads,mode):
        if mode=="mem":
#            print(int(mem), wildcards)
            return int(mem)
        if mode=="partition":
#            print(partition,wildcards)
            return partition
        if mode=="threads":
            return int(threads)

    # change with each attempt
    # Where input.size//1000000 is used convert the cumulative size of input files in bytes to mb, and the tailing 2 could be any arbitrary number based on the specifics of your shell/script requirements.
    mem = max((input.size//1000000) * attempt * mult, attempt*min_size* mult) # this is mb

    # handle case where we are not on a cluster, no partition is defined
    if SLURM_PARTITIONS[0][0]=="":
        partition = ""
        return return_result(mem,partition,threads,mode)

    # choose the best partition, priority to memory
    # HIGH_MEM_PARTITIONS = [["name","min_memory","max_memory","min_threads","max_threads"]]

    #### check memory, get all partition giving mem/1000 (Gb)
    mem_selection = [part for part in SLURM_PARTITIONS if part[2]>=mem]

    # if empty, get the partition with most mem
    if not mem_selection:
        mem_selection = [max(SLURM_PARTITIONS,key=lambda x:x[2])]

    #### check threads, get all partitions giving at least threads
    threads_selection = [part for part in mem_selection if part[4]>=threads]

    # if empty, get partition with most threads
    if not threads_selection:
        threads_selection = [max(mem_selection,key=lambda x:x[4])]

    #### if more than 1 partition remain, 
    # choose the one with, in this order 
    # least min_mem, min_threads, max_mem, max_threads
    final_selection = min(threads_selection,key=lambda x:[x[1],x[3],x[2],x[4]])

    #### decide on mem/threads to use
    partition, min_memory, max_memory, min_threads ,max_threads = final_selection
    mem_final = min(max(mem, min_memory), max_memory)
    threads_final = min(max(threads, min_threads), max_threads)

    # output
    return return_result(mem_final,partition,threads_final,mode)






    
