#!/usr/bin/env python3
import os
from os.path import dirname
import argparse

# get strong env path..... for reason
def which(pgm):
    # https://stackoverflow.com/questions/9877462/is-there-a-python-equivalent-to-the-which-command
    path=os.getenv('PATH')
    for p in path.split(os.path.pathsep):
        p=os.path.join(p,pgm)
        if os.path.exists(p) and os.access(p,os.X_OK):
            return p


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path_gtdb", help="path to gtdb database")
    args = parser.parse_args()

    GTDB = args.path_gtdb

    # get strong env path
    path = which("gtdbtk")
    env_path = path.split("/bin/gtdbtk")[0]
    # just check it is not empty
    if len(os.listdir(GTDB))<=1:
        # create a copy of download script with new dl directory
        dl_script = which("download-db.sh")
        new_dl_script = dirname(dl_script)+"/download-db2.sh"
        os.system("cp %s %s"%(dl_script,new_dl_script))
        os.system("sed -i 's=${{GTDBTK_DATA_PATH}}=%s=g' %s"%(GTDB,new_dl_script))
        os.system("download-db2.sh &>{log}")
    # overwrite default database location
    export_path = env_path+"/etc/conda/activate.d/gtdbtk.sh"
    with open(export_path,"w") as handle:
        handle.write("export GTDBTK_DATA_PATH=%s/\n"%GTDB)
