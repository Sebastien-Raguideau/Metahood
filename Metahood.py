#!/usr/bin/env python3
from os.path import abspath, realpath, dirname, basename, exists
from scripts.common import fill_default_values, cd 
from subprocess import PIPE,Popen
from psutil import virtual_memory
import subprocess
import argparse
import shutil
import yaml
import time
import sys
import os



parser = argparse.ArgumentParser(description="MetaHood - pipeline for assembly and binning of metagenomic samples")
parser.add_argument("step", nargs='?', default='all', type=str, choices={"all", "sample_qc"}, help="Pipeline step to run options: all or sample_qc, default is all")
parser.add_argument("config", type=str, help="config_file.yaml to use")
parser.add_argument("--cores", "-c", type=int, default=1, help="Number of threads")
parser.add_argument("--verbose", "-v", action="store_true", help="Increase verbosity level")
parser.add_argument("--dryrun", "-n", action="store_true", help="Show tasks, do not execute them")
parser.add_argument("--unlock", "-u", action="store_true", help="Unlock the directory")
parser.add_argument("--touch", "-t", action="store_true", help="Touch all files, to reset timestamp and stop unwanted/uncalled snakemake reruns")
parser.add_argument("--dag", "-d", help="file where you want the dag to be stored")
parser.add_argument('-s', nargs=argparse.REMAINDER,help="Pass additional argument directly to snakemake")
args = parser.parse_args()

# get config file
CONFIG_FILE = abspath(realpath(args.config))
config = yaml.full_load(open(CONFIG_FILE))

# get exec directory
METAHOOD_DIR = dirname(abspath(realpath(sys.argv[0])))

# execution directory
EXEC_DIR=abspath(realpath(config["execution_directory"]))
os.system("mkdir -p %s"%EXEC_DIR)

# ------- set max memory used, in Go ---------------
Mem_tot=virtual_memory().total
Percent_mem=config["Percent_memory"]
MEMG=str(int((Percent_mem*Mem_tot)/10**9))

# ------- base parameters used to call snakemake -----------
base_params = ["snakemake", "--directory", EXEC_DIR, "--cores", str(args.cores), "--config", "LOCAL_DIR=%s"%METAHOOD_DIR,"CONFIG_PATH=%s"%CONFIG_FILE,"EXEC_DIR=%s"%EXEC_DIR,"--configfile="+CONFIG_FILE,"--resources",'memG='+MEMG, "--latency-wait", "120","-k","--use-conda"]

# ------- additional parameters -----------
if args.verbose:
    base_params.extend(["-p", "-r", "--verbose"]) 
if args.dryrun:
    base_params.extend(["--dryrun"])
if args.unlock:
    base_params.extend(["--unlock"])
if args.touch:
    base_params.extend(["-t"])    
if args.dag:
    base_params.extend(["--rulegraph"])
if args.s :
    base_params.extend(args.s)


# ------- call snakemake from  METAHOOD_DIR -----------
with cd(METAHOOD_DIR):
    def call_snake(extra_params=[]):
        call_snake.nb+=1
        if args.dag:
            p1=Popen(base_params + extra_params, stdout=PIPE, stderr=sys.stderr)
            p2=Popen(["dot","-Tpng"],stdin=p1.stdout, stdout=PIPE, stderr=sys.stderr)
            with open(args.dag.replace(".png",str(call_snake.nb)+".png"),"bw") as f :
                f.write(p2.communicate()[0])
        else :
            subprocess.check_call(base_params + extra_params, stdout=sys.stdout, stderr=sys.stderr)
    call_snake.nb=0
    if args.step == 'sample_qc':
        call_snake(["--snakefile", "sample_qc.snake"])
    if args.step == 'all':
        #launch master snake
        call_snake(["--snakefile", "Master.snake"])
    # setup data folder
    # call_snake(["--snakefile", "Setup_samples.snake"])
    # launch Maganalysis
    # call_snake(["--snakefile", "Maganalysis.snake"])    
    # launch desman
    # call_snake(["--snakefile", "Desman.snake"])


