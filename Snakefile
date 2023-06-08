# run on login node / compute server
#snakemake -p -j 3

# run on compute node of HPC
#snakemake -p --cluster-sync "srun --cpus-per-task=6 --mem-per-cpu=1875 --time=7:59:00" --jobs 100

# python module import
import glob
import getpass

# start snakemaking

# get parameters via config file
configfile: "config.yaml"

# directories
RESDIR = config["BASEDIR"]+"/"+config["VIRFAM"]+"/"+config["PROJECTID"]+"/results"
LOGDIR = config["BASEDIR"]+"/"+config["VIRFAM"]+"/"+config["PROJECTID"]+"/logs"
TABDIR = config["BASEDIR"]+"/"+config["VIRFAM"]+"/"+config["PROJECTID"]+"/hittables"
if not os.path.exists( config["FASTQDIR"] ):
	os.system( "mkdir "+config["FASTQDIR"] )
if not os.path.exists( TABDIR ):
	os.system( "mkdir "+TABDIR )

# check databases
import subprocess

for db in (config["DBFILTER"], config["DBREFSEQ"], config["DBVIRAL"]):
	command = ['blastdbcmd', '-db', db, '-info']
	p = subprocess.run(command, stdout=subprocess.PIPE)

if not os.path.isfile( config["ACCSVIRAL"] ):
	raise AssertionError(f'\nA file with viral protein accessions is missing or not a file. '
                         f'Please, check the file: {config["ACCSVIRAL"]}\n')

# process ID and username
PID = os.getpid()
USR = getpass.getuser()

# temporary directory
TDIR= "/tmp/"+str(USR)+"-"+str(PID)

# sample IDs
SAMPLES = []
if config["ACCLIST"] == "":
	FILES   = [ os.path.basename(x)    for x in glob.glob( config["FASTQDIR"]+"/*.fastq.gz" ) ]
	FILES2  = [ os.path.splitext(x)[0] for x in FILES ]
	SAMPLES = [ os.path.splitext(x)[0] for x in FILES2 ]
else:
	os.system( "mkdir "+TDIR )
	os.system( "wd=`pwd`; cd "+TDIR+"; awk \'{printf \"%s\\n\", $1>$1\".txt\"}\' "+config["ACCLIST"]+"; cd $wd" )
	IDFILES = [ os.path.basename(x)    for x in glob.glob( TDIR+"/*.txt" ) ]
	for idf in IDFILES:
		if not os.path.exists( config["FASTQDIR"]+"/"+idf ):
			os.system( "cp "+TDIR+"/"+idf+" "+config["FASTQDIR"]+"/" )
	os.system( "rm -rf "+TDIR )	
	FILES   = [ os.path.basename(x)    for x in glob.glob( config["FASTQDIR"]+"/*.txt" ) ]
	SAMPLES = [ os.path.splitext(x)[0] for x in FILES ]

# for debugging	
#print(RESDIR)
#print( ",".join(SAMPLES) )
#print( ",".join(SAMPLES) )


# start calculating

# run everything
rule all:
 message:
  "Running everything"
 input:
  expand( config["FASTQDIR"]+"/{sample}.fastq.gz",                                   sample=SAMPLES ),
  expand( RESDIR+"/{sample}/virushunter/contigs.singlets.fas.gz",                    sample=SAMPLES ),
  expand( RESDIR+"/{sample}/virusgatherer/genseedhmm-"+config["ASSEMBLER"]+".fasta", sample=SAMPLES )


# virushunter search
rule hunter:
 message:
  "Executing virushunter search"
 input:
  config["FASTQDIR"]+"/{sample}.fastq.gz"
 output:
  RESDIR+"/{sample}/virushunter/contigs.singlets.fas.gz"
 log:
  err   = LOGDIR+"/virushunter/{sample}_virushunter.err"
 params:
  dir1  = config["WFLOWDIR"],
  dir2  = config["BASEDIR"],
  vfam  = config["VIRFAM"],
  cpus  = config["THREADS"],
  sid   = "{sample}",
  pid   = config["PROJECTID"],
  db1   = config["DBFILTER"],
  db2   = config["DBREFSEQ"],
  db3   = config["DBVIRAL"],
  db4   = config["ACCSVIRAL"],
  flag1 = config["MAPHG38"],
  flag2 = config["DEBUGMODE"]
 shell:
  "{params.dir1}/1_scripts/virushunterTWC.pl {params.vfam} {input} {params.pid} {params.sid} {params.dir2} {params.cpus} {params.db1} {params.db2} {params.db3} {params.db4} {params.dir1} {params.flag1} {params.flag2} 2> {log.err}"


# virushunter hit table
rule hunter_hittab:
 message:
  "Collecting virushunter results"
 input:
  expand( RESDIR+"/{sample}/virushunter/contigs.singlets.fas.gz", sample = SAMPLES )
 output:
  TABDIR+"/virushunter.tsv"
 params:
  dir1 = config["WFLOWDIR"],
  dir2 = config["BASEDIR"],
  vfam = config["VIRFAM"],
  pid  = config["PROJECTID"],
  flag = config["ISSRA"],
 shell:
  "{params.dir1}/1_scripts/virushunterTWC_hittable.pl {params.vfam} {params.pid} {params.dir2} {params.flag} > {output}"


# download data if not present
rule SRA_download: 
 message:
  "Downloading SRA data"
 input:
  config["FASTQDIR"]+"/{sample}.txt"
 output:
  config["FASTQDIR"]+"/{sample}.fastq.gz"
 params:
  dir1 = config["WFLOWDIR"],
  dir2 = config["FASTQDIR"]
 shell:
  "{params.dir1}/1_scripts/downloadFromSRA.pl {input} {params.dir2}"


# virusgatherer assembly
rule gatherer:
 message:
  "Executing virusgatherer assembly using "+config["ASSEMBLER"]
 input:
  fastq = config["FASTQDIR"]+"/{sample}.fastq.gz",
  seed  = RESDIR+"/{sample}/virushunter/contigs.singlets.fas.gz"
 output:
  RESDIR+"/{sample}/virusgatherer/genseedhmm-"+config["ASSEMBLER"]+".fasta"
 log:
  err   = LOGDIR+"/virusgatherer/{sample}_virusgatherer.err"
 params:
  dir1  = config["WFLOWDIR"],
  dir2  = config["BASEDIR"],
  vfam  = config["VIRFAM"],
  cpus  = config["THREADS"],
  sid   = "{sample}",
  pid   = config["PROJECTID"],
  ass   = config["ASSEMBLER"],
  db1   = config["DBREFSEQ"],
  db2   = config["ACCSVIRAL"],
  flag1 = config["DEBUGMODE"]
 shell:
  "{params.dir1}/1_scripts/virusgathererTWC.pl {params.vfam} {input.fastq} {params.pid} {params.sid} {input.seed} {params.dir2} {params.cpus} {params.ass} {params.db1} {params.db2} {params.dir1} {params.flag1} 2> {log.err}"



