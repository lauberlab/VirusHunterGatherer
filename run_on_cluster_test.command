
snakemake -n -p -j 3 --latency-wait 60 --cluster-config cluster.yaml --cluster "bsub -R 'rusage[mem={cluster.mem}] affinity[core({cluster.cores})]' -W {cluster.time} -o test.out -e test.err {cluster.mail}"

