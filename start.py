#!/usr/bin/env python3
import argparse
import subprocess
from subprocess import PIPE,Popen
from psutil import virtual_memory
import sys
import os
import os.path
import shutil
import yaml
import time

from scripts.common import fill_default_values


parser = argparse.ArgumentParser(description="MetaHood - pipeline for assembly and binning of metagenomic samples")
parser.add_argument("--threads", "-t", type=int, default=1, help="Number of threads")
parser.add_argument("dir", type=str, help="Output directory")
parser.add_argument("--config", "-c", type=str, default="", help="config_file.yaml to use")
parser.add_argument("--verbose", "-v", action="store_true", help="Increase verbosity level")
parser.add_argument("--dryrun", action="store_true", help="Show tasks, do not execute them")
parser.add_argument("--unlock", "-u", action="store_true", help="Unlock the directory")
parser.add_argument("--dag", "-d", help="file where you want the dag to be stored")
parser.add_argument('-s', nargs=argparse.REMAINDER,help="Pass additional argument directly to snakemake")
args = parser.parse_args()

exec_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
LOCAL_DIR = os.path.realpath(exec_dir)
print(LOCAL_DIR)

# ------- set max memory used, in Go ---------------
CONFIG_FILE = os.path.abspath(args.config)
Mem_tot=virtual_memory().total
Percent_mem=yaml.full_load(open(CONFIG_FILE))["Percent_memory"]
MEMG=str(int((Percent_mem*Mem_tot)/10**9))

# ------- base parameters used to call snakemake -----------
base_params = ["snakemake", "--directory", os.path.realpath(args.dir), "--cores", str(args.threads), "--config", "LOCAL_DIR" + "=" + LOCAL_DIR,"--configfile="+CONFIG_FILE,"--resources",'memG='+MEMG, "--latency-wait", "12","--use-conda"]

if args.verbose:
    base_params.extend(["-p", "-r", "--verbose"]) 
if args.dryrun:
    base_params.extend(["--dryrun"])
if args.unlock:
    base_params.extend(["--unlock"])
if not os.path.exists(args.dir):
    os.makedirs(args.dir)
if args.dag:
    base_params.extend(["--rulegraph"])
if args.s :
    base_params.extend(args.s)


# ------- call snakemake from  exec_dir -----------

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


with cd(exec_dir):
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
    # setup data folder
    call_snake(["--snakefile", "Setup_samples.snake"])
    #launch master snake
    call_snake(["--snakefile", "Master.snake"])
    # launch desman
    call_snake(["--snakefile", "Desman.snake"])


