from __future__ import print_function
from __future__ import division

try:
    from future_builtins import zip
except:
    pass

from collections import defaultdict
from subprocess import Popen, PIPE
from os.path import basename
import os
import os.path
import re


default_values = {
    "binning":{"concoct":{"contig_size" : 1000,"execution" : 1,"max_bin_nb" : 2000},"metabat2":{"execution" : 1,"contig_size":1500}},
    "mag":["native"],
    "threads":8,
    "assembly":    {"assembler": "megahit","groups": {},"parameters":"" },
    "annotation": {'diamond':{},"cat_db":"","kraken_db":"","kofamscan":{"profiles":"","ko_list":""}},
    "graph":{"List_graphs":{}},
    "filtering":"",
    "Percent_memory":0.5,
    "task_memory":200,
    "maganalysis":0,
    "desman":{"execution":0, "nb_haplotypes": 10,"nb_repeat": 5,"min_cov": 1,"scripts":""}
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
        default_values["scg_data"]= os.path.join(local_dir, "scg_data")
        default_values["conda_env"]= os.path.join(local_dir, "conda_envs")
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

def replace_extensions(sample):
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
