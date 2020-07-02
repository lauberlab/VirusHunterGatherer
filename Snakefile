# run on login node / compute server
#snakemake -p -j 3

# run on compute node of HPC
#snakemake -p --cluster-sync "srun --cpus-per-task=6 --mem-per-cpu=1875 --time=7:59:00" --jobs 100

# python module import
import glob

# start snakemaking

# get parameters via config file
configfile: "config.yaml"

# process input
FILES   = [ os.path.basename(x)    for x in glob.glob( config["FASTQDIR"]+"/*.fastq.gz" ) ]
FILES2  = [ os.path.splitext(x)[0] for x in FILES ]
SAMPLES = [ os.path.splitext(x)[0] for x in FILES2 ]
RESDIR  = config["BASEDIR"]+"/"+config["VIRFAM"]+"/"+config["PROJECTID"]+"/results"

# start calculating
rule all:
 message:
  "Running everything"
 input:
  expand( RESDIR+"/{sample}/final.hits.tsv", sample=SAMPLES )

rule hunter:
 message:
  "Executing virushunter search"
 input:
  fastq=config["FASTQDIR"]+"/{sample}.fastq.gz",
 output:
  RESDIR+"/{sample}/final.hits.tsv"
 params:
  dir1=config["WFLOWDIR"],
  dir2=config["BASEDIR"],
  vfam=config["VIRFAM"],
  cpus=config["THREADS"],
  sid="{sample}",
  pid=config["PROJECTID"],
  db1=config["DBFILTER"],
  db2=config["DBREFSEQ"],
  db3=config["DBVIRAL"],
  db4=config["ACCSVIRAL"],
  flag1=config["MAPHG38"],
  flag2=config["DEBUGMODE"]
 shell:
  "{params.dir1}/1_scripts/virushunterTWC.pl {params.vfam} {input.fastq} {params.pid} {params.sid} {params.dir2} {params.cpus} {params.db1} {params.db2} {params.db3} {params.db4} {params.dir1} {params.flag1} {params.flag2}"

