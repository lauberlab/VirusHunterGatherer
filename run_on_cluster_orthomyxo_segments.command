#!/bin/bash

snakemake -n -k -p -j 20 --latency-wait 60 --cluster-config cluster.yaml --cluster "bsub -R 'rusage[mem={cluster.mem}] affinity[core({cluster.cores})]' -W {cluster.time} -o orthomyxo.log1 -e orthomyxo.log1 {cluster.mail}"

