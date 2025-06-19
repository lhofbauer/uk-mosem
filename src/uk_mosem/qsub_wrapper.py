"""
Qsub wrapper script for snakemake


Copyright (C) 2025 Leonhard Hofbauer, licensed under a MIT license

"""

import os
import sys

from snakemake.utils import read_job_properties

js = sys.argv[1]
jp = read_job_properties(js)

# format params
threads = jp["threads"]
runtime = '{:02d}:{:02d}:00'.format(*divmod(jp["resources"]["runtime"], 60)) 
mem = jp["resources"]["mem_mb"]/1024
# access property defined in the cluster configuration file (Snakemake >=3.6.0)
# job_properties["cluster"]["time"]

os.system("qsub -l h_rt={runtime} -l mem={mem}G  -pe smp {threads} -m be {script}".format(runtime=runtime,
                                                                                          mem=mem,
                                                                                          threads=threads,
                                                                                          script=js))
