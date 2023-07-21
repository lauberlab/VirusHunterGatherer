
snakemake -n -k -p -j 20 --latency-wait 60 --cluster-config cluster.yaml --cluster "bsub -R 'rusage[mem={cluster.mem}] affinity[core({cluster.cores})]' -W {cluster.time} -o flavi.out -e flavi.err {cluster.mail}"

