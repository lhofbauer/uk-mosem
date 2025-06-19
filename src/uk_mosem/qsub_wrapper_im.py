"""
Qsub wrapper script for snakemake (including for immediate submits)


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license

"""

import os
import sys
import random
from snakemake.utils import read_job_properties

js = sys.argv[-1]
dependencies = set(sys.argv[1:-1])

jp = read_job_properties(js)

# format params
threads = jp["threads"]
runtime = '{:02d}:{:02d}:00'.format(*divmod(jp["resources"]["runtime"], 60)) 
mem = max(1,round(jp["resources"]["mem_mb"]/1024))
# access property defined in the cluster configuration file (Snakemake >=3.6.0)
# job_properties["cluster"]["time"]

d="./logs/"
ri = random.randint(0, 1000)
subline = ["qsub -o {outlog} -e {errlog}".format(outlog=d+jp["rule"]+str(ri)+"_out.txt",
                                                 errlog=d+jp["rule"]+str(ri)+"_err.txt")]

if dependencies:
    #subline.append("--dependency")
    # only keep numbers in dependencies list
    #dependencies = [ x for x in dependencies if x.isdigit()]
    subline.append("-hold_jid " + ",".join(dependencies))

if True:#jp["rule"] == "all":
    args="-l h_rt={runtime} -l mem={mem}G  -pe smp {threads} -m e {script}".format(runtime=runtime,
                                                                                    mem=mem,
                                                                                    threads=threads,
                                                                                    script=js)
    # -m be
else:
    args="-l h_rt={runtime} -l mem={mem}G  -pe smp {threads} {script}".format(runtime=runtime,
                                                                                    mem=mem,
                                                                                    threads=threads,
                                                                                    script=js)
subline.append(args)


# adjust as needed
subline.append(" | cut -c 10-16")


os.system(" ".join(subline))

