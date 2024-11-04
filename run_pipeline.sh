#!/bin/bash

snakemake \
--snakefile Snakefilefilter_keepRPMDPM.smk \
--use-conda \
--rerun-incomplete \
-j 64 \
--cluster-config cluster_filter.yaml \
--configfile config.yaml \
--cluster "sbatch -c {cluster.cpus} \
-t {cluster.time} -N {cluster.nodes} \
--mem {cluster.mem} \
--output {cluster.output} \
--error {cluster.error}"
