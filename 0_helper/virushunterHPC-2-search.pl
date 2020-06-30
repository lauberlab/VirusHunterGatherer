#!/usr/bin/perl

#########################################################
# search for viral sequences in SRA
#  - retrieve SRA IDs from NCBI taxonomy info
#  - retrieve SRA reads through fastq-dump (SRA toolkit)
#  - tblastn/hmmsearch against specified viral query
#
# author: Chris Lauber
#########################################################

# ------------ #
# load modules
# ------------ #
use lib "/home/lauber/lib/perl/";
use warnings;
use strict;
use LWP::Simple qw/get/;
use POSIX qw/floor/;
use Cwd;
use Data::Dumper;
use Parallel::ForkManager;

use lib "/home/lauber/lauber-2015-virusHunter";
use      virusHGutils;


# --------- #
# parameter
# --------- #
my $resdir        = "/projects/p_sra";
my $tmpdir1       = "/tmp";
my $tmpdir2       = sprintf "%s/virushunter",             virusHGutils::get_workspace_path( "lustre" );
my $tmpdir3       = sprintf "%s/lauber-2015-virusHunter", virusHGutils::get_workspace_path( "scratch" );
my $tmpdir;
my $logtmpdir     = sprintf "%s/virushunter",             virusHGutils::get_workspace_path( "lustre" );
my $ncbirootASCP  = "130.14.29.30";
my $ncbirootFTP   = "130.14.250.7";
my $asperaroot    = "~/.aspera/connect";
my $sampleID      = "";
my $sraSize       = 0;
my $blastE        = 1;
my $blastEsure    = 1e-04;
my $hmmerE        = 10;
my $hmmerEsure    = 0.01;
my $identTH       = 50;
my $identTHn      = 2;
my $motifF        = "";
my $maxHitN       = 10000;
my $readHitLenTH  = 20;
my $hitsPerReadTH = 2;
my $hitsTotalTH   = 1;
my $hitsMotifsTH  = 1;
my $hitsTotalHmm  = 2;
my $fimoTHq       = 0.001;
my $threads       = 1;
my $debugMode     = 0;
my $doProfile     = 0;
my $doFilter      = 0;
my $doRefSeqBlast = 1;
my $refseqDB      = sprintf "%s/db/RefSeq/refseq_protein",        virusHGutils::get_workspace_path( "scratch" );
my $viralDB       = sprintf "%s/db/RefSeq/viral_genomic",         virusHGutils::get_workspace_path( "scratch" );
my $refseqNegGIs  = sprintf "%s/db/RefSeq/viral_protein.gi_list", virusHGutils::get_workspace_path( "scratch" );
my $filterDB      = sprintf "%s/db/RefSeq/filter",                virusHGutils::get_workspace_path( "scratch" );
#my $filterDB      = sprintf "%s/db/Kroepfchen/Kroepfchen_050_entkarpft_reformatted_filtered", virusHGutils::get_workspace_path( "scratch" );
my $transcriptsHg38 = sprintf "%s/db/genomes/human/hg38/hg38-cds-index",           virusHGutils::get_workspace_path( "scratch" );
my $genemappingHg38 = sprintf "%s/db/genomes/human/hg38/hg38-cds-geneMapping.tsv", virusHGutils::get_workspace_path( "scratch" );
my $filterHg38    = 0;
my $mapHg38       = 0;
my $gi2taxExec    = "/home/lauber/lauber-2015-virusHunter/virusutils-queryTax.pl";
my $DBblastEh     = 1e-04;	# cut-off for blast against host sequences
my $DBblastEv     = 1; # 1e-02;	# cut-off for blast against virus sequences
my $RefSeqBlastMx = 5;
my $cutQth        = 10;	# threshold for cutting low-quality bases in cutadapt
my $assembler     = "cap3";
my $cap3overhang  = 75; # in %
my $cap3overlap   = 20; # in bases
my $cap3ident     = 85; # in %
my $filterFile    = "";
my $filterE       = 1e-05;
my $filterIdent   = 99;
my $sizeInflation = 20; # factor for estimating total required size from size of initial .sra file
my $fastq_input   = 0;
my $threadsMax    = 32; # don't use more threads than this (to run multiple jobs on very large nodes)

# --------------- #
# usage and input
# --------------- #
if ( $#ARGV < 3 ){
die("
usage: virushunterHPC-2-search.pl <family> <query_fasta_file> <project_id> <sra_run_id> [options]\n
options:
\t-m=<string>\tsearch for motifs given in file <m> in MEME format
\t-e=<float>\tE-value threshold during search
\t\t\t(default: $blastE for tblastn, $hmmerE for hmmsearch)
\t-t=<integer>\tnumber of threads used for tblastn (default: $threads)
\t-r=<integer>\trequire at least <r> read hits (default: $hitsTotalTH)
\t-h=<integer>\tconsider as hit if at least one read hits at least
\t\t\t<h> queries (default: $hitsPerReadTH)
\t-i=<float>\tconsider as hit if <j> (see below) read hits have
\t\t\tidentity of at least <i> (default: $identTH)
\t-j=<integer>\tconsider as hit if <j> read hits have identity of 
\t\t\tat least <i> (see above) (default: $identTHn)
\t-l=<integer>\tminimum length of a read hit [aa] (default: $readHitLenTH)
\t-n=<integer>\tonly report best <n> read hits (default: $maxHitN)
\t-f=<string>\trun additional filtering step against known viral
\t\t\tgenomes provided in file <f>
\t-p\t\tdo profile-based HMMer search (default: Blast search)
\t-rsb\t\tdo not run counter-blast search against refseq_protein
\t-fastq\t\tprovided are not SRA identifiers for download but names of
\t\t\tlocal fastq.gz files
\t-hg38f\t\tfilter against hg38 in filtering step 1
\t-hg38m\t\tmap reads against hg38 genes (only meaningful for transcriptomes)
\t-d\t\tdebug mode; temporary files will be kept
\n");
}

# -------------------- #
# log program call
# -------------------- #
my $virushuntercall = "virushunterHPC-2-search.pl ".join(" ",@ARGV);
print "$virushuntercall\n\n";

# -------------------- #
# get input parameters
# -------------------- #
my $family = shift;
my $queryF = shift;
my $projID = shift;
my $sraid  = shift;
while ( $#ARGV > -1 ){
	my $option = shift;
	if ( $option =~ /-m=(.*)/  ){ $motifF         = $1; }
	if ( $option =~ /-e=(.*)/  ){ $blastE         = $1; $hmmerE = $1; }
	if ( $option =~ /-t=(\d+)/ ){ $threads        = $1; }
	if ( $option =~ /-r=(\d+)/ ){ $hitsTotalTH    = $1; }
	if ( $option =~ /-h=(\d+)/ ){ $hitsPerReadTH  = $1; }
	if ( $option =~ /-i=(.*)/  ){ $identTH        = $1; }
	if ( $option =~ /-j=(\d+)/ ){ $identTHn       = $1; }
	if ( $option =~ /-l=(\d+)/ ){ $readHitLenTH   = $1; }
	if ( $option =~ /-n=(\d+)/ ){ $maxHitN        = $1; }
	if ( $option =~ /-f=(.+)/  ){ $filterFile     = $1; $doFilter = 1; }
	if ( $option eq "-p"       ){ $doProfile      =  1; }
	if ( $option eq "-rsb"     ){ $doRefSeqBlast  =  0; }
	if ( $option eq "-fastq"   ){ $fastq_input    =  1; }
	if ( $option eq "-hg38f"   ){ $filterHg38     =  1; }
	if ( $option eq "-hg38m"   ){ $mapHg38        =  1; }
	if ( $option eq "-d"       ){ $debugMode      =  1; }
}


# -------------- #
# initialization
# -------------- #
my $transeqbin = "transeq";
my ($stime, $rtimeD, $rtimeP, $rtimeS, $rtimeF, $rtimeM) = (0,0,0,0,0,0);
$stime = time();
init();

# ------------ #
# prepare data
# ------------ #
preprocess();

# ------------------------------------- #
# optional filter against known viruses
# ------------------------------------- #
my %filterHits = ();
filter() if ( $doFilter == 1 );

# -------------------------------- #
# run Blast or HMMer against reads
# -------------------------------- #
if ( $doProfile == 1 ){
	searchByHMMer();
}else{
	searchByBlast();
}

# ------------------------------------ #
# run Blast against reference proteins
# ------------------------------------ #
if ( $doRefSeqBlast == 1 ){
	searchRefSeq();
}

# ----------------------------------- #
# map against reference transcriptome
# ----------------------------------- #
if ( $mapHg38 == 1 ){
	mapRefGenes();
}

# ------------------- #
# compress some files
# ------------------- #
compress();

# --------------------------- #
# total runtime and data size
# --------------------------- #
reportResources();

# ----- #
# clean
# ----- #
cleanup();

# --------- #
# functions
# --------- #
# init project (if started the first time)
sub init {
	# set result and temporary directories
	$resdir    .= "/$family/results";
	$logtmpdir .= "/$family/$projID/log";
	$tmpdir     = "$tmpdir3/$family";

	# set databases
	$viralDB   .= "_$family"  if ( -e $viralDB."_$family.nhr"  or  -e $viralDB."_$family.nal" );
	$filterDB  .= "_hg38"     if $filterHg38 == 1;

	# was this SRA run already processed before?
	open (OLD, "<$resdir/$projID/log/completed.txt" );
	while ( my $line = <OLD> ){
		chomp($line);
		next if ( $line eq "" );
		exit() if $line eq $sraid;
	}
	close(OLD);

	# delete all own files on local SSD, as a precaution
	`find $tmpdir1/ -user lauber -exec rm -fr {} \\; 2> /dev/null`;

	# determine size requirements for this run
	my $srasize = 0;
	my $srafile_name = "$sraid.sra";
	   $srafile_name = "$sraid.fastq.gz"  if ( $fastq_input );
	$srasize = `du $tmpdir/$srafile_name | awk '{print \$1};'`;  chomp( $srasize );

	# reset temporary data directory depending on size requirements
	my $diskfree1 = `df $tmpdir1 | tail -n 1 | awk '{print \$4};'`;  chomp( $diskfree1 ); # free space on  local SSD of compute node
	my $diskfree2 = `df $tmpdir2 | tail -n 1 | awk '{print \$3};'`;  chomp( $diskfree2 ); # free space on global SSD on /scratch
	$diskfree2 -= 5 * 1024 * 1024 * 1024; # always keep 5 Terabyte free on global SSD
	$tmpdir = "$tmpdir2/$family" if ( $srasize * $sizeInflation < $diskfree2 );
	$tmpdir = "$tmpdir1/$family" if ( $srasize * $sizeInflation < $diskfree1 );
	if ( $tmpdir ne "$tmpdir3/$family" ){
		mkdir( "$tmpdir" )         if ! -e "$tmpdir";
		mkdir( "$tmpdir/$projID" ) if ! -e "$tmpdir/$projID";
		`cp     $tmpdir3/$family/$srafile_name $tmpdir/`;
		`rm -rf $tmpdir3/$family/$srafile_name`;
	}

	# report hostname and data directory
	my $hostname = `hostname`; chomp( $hostname );
	printf "[virushunter] \t$sraid: running on %s\n", $hostname;
	printf "[virushunter] \t$sraid: using temporary data directory %s\n\n", $tmpdir;

	# verify that transeq binary in path
	my $whichtranseq = `which transeq 2>/dev/null`;
	if ( $whichtranseq eq "" ){
		$transeqbin = "/home/lauber/lib/EMBOSS-6.5.7/emboss/transeq";
	}

	# limit threads to 32
	$threads = $threadsMax       if ( $threads >= ($threadsMax*2) );
	$threads = int($threads / 2) if ( $threads <= ($threadsMax*2) and $threads > $threadsMax );
	
} # end init

# clean temporary data
sub cleanup {
	printf "\n[virushunter] \t$sraid: cleaning temporary data and finishing up\n";
	# log this SRA run as done
	open( PROGRESS, ">>$tmpdir3/$family/virushunterHPC-$projID-progress.txt" ) or die("$0: could not open file: $!\n");
	flock(PROGRESS,2) or die("$0: could not lock file: $!\n");
	print PROGRESS $sraid."\n";
	close(PROGRESS) or die("$0: could not close file: $!\n");
	
	# log global statistics
	open(  STATS, ">>$tmpdir3/$family/virushunterHPC-$projID-search.txt" ) or die("$0: could not open file: $!\n");
	flock( STATS, 2 ) or die("$0: could not lock file: $!\n");
	printf STATS "%s\t%s\t%.4f\t%d\t%d\n", $sraid, $projID, $sraSize, $rtimeS+$rtimeF, $threads;
	close( STATS ) or die("$0: could not close file: $!\n");

	# move result file to /scratch if not already there
	if ( $tmpdir ne "$tmpdir3/$family" ){
		`cp -r  $tmpdir/$projID/$sraid $tmpdir3/$family/$projID/`;
		`rm -rf $tmpdir/$projID/$sraid`;
		`rm -rf $tmpdir`  if $tmpdir eq "$tmpdir1/$family";
	}

	# delete temporary files
	if ( $debugMode == 0  or  $tmpdir ne "$tmpdir3/$family" ){
		#`rm $logtmpdir/$sraid-data.txt`;
		`rm -rf $tmpdir/$sraid*`;
	}
	if ( $debugMode == 0 ){
		`rm -rf $tmpdir3/$family/$sraid*.slurm`;
	}

	# log this runs as completed
	open( LOG, ">>$tmpdir3/$family/virushunterHPC-$projID-completed.txt" ) or die("$0: could not open file: $!\n");
	flock(LOG, 2) or die("$0: could not lock file: $!\n");
	print LOG $sraid."\n";
	close(LOG) or die("$0: could not close file: $!\n");
	
	# THE END
	printf "[virushunter] \t$sraid: done\n\n";
} # end cleanup

# preprocess data for Blast or HMMer search
sub preprocess {
	# verify that SRA dataset was successfully downloaded
	exit() if ( $fastq_input == 0 and ! -e "$tmpdir/$sraid.sra" );
	# transform to Fasta format
	if ( ! -e "$tmpdir/$sraid.fasta" ){
		printf "[virushunter] \t$sraid: transforming to Fasta format\n";
		#`fastq-dump -O $tmpdir --fasta 0 -B $tmpdir/$sraid.sra`;
		#`fastq-dump -O $tmpdir --fasta 0 -B --split-spot --skip-technical --readids $tmpdir/$sraid.sra`;
		# uncomment for fastq file generation
		#`fastq-dump -O $tmpdir -B --split-spot --skip-technical --readids $tmpdir/$sraid.sra`;
		#`gzip $tmpdir/$sraid.fastq`;
		#`mv   $tmpdir/$sraid.fastq.gz $tmpdir3/$family/`;
		if ( $fastq_input == 0 ){
			`fastq-dump -O $tmpdir -B --split-spot --skip-technical --readids --gzip --clip $tmpdir/$sraid.sra`;
			#`fasterq-dump $tmpdir/$sraid.sra -O $tmpdir -t $tmpdir --split-spot --skip-technical --threads $threads`;
		}
		my $read1id = `zcat $tmpdir/$sraid.fastq.gz | head -n 1`;
		my $read2id = `zcat $tmpdir/$sraid.fastq.gz | head -n 5 | tail -n 1`;
		chomp($read1id);
		chomp($read2id);
		$read1id =~ /(.*\.\d+)\.\d .*/;  my $r1id = $1;
		$read2id =~ /(.*\.\d+)\.\d .*/;  my $r2id = $1;
		my $fastp_w   = $threads;
		   $fastp_w   = 16  if $threads > 16;
		my $fastp_cmd = "fastp -i $tmpdir/$sraid.fastq.gz -w $fastp_w -j $tmpdir/$sraid.fastp.json -h $tmpdir/$sraid.fastp.html";
		if ( $r1id eq $r2id ){
			$fastp_cmd .= " --interleaved_in --stdout > $tmpdir/$sraid.trim.fastq";
		}else{
			$fastp_cmd .= " -o $tmpdir/$sraid.trim.fastq";
		}
		`$fastp_cmd`;
		`seqtk seq -A $tmpdir/$sraid.trim.fastq > $tmpdir/$sraid.fasta`;
		`rm $tmpdir/$sraid.fastq.gz`;
		`rm $tmpdir/$sraid.trim.fastq`;
		#exit;
	}
	# create Blast database if needed
	if ( $doProfile == 0  or  $doFilter == 1 ){
		if ( ! -e "$tmpdir/$sraid.fasta.nhr" and ! -e "$tmpdir/$sraid.fasta.00.nhr" ){
			printf "[virushunter] \t$sraid: creating blast DB\n";
			`makeblastdb -in $tmpdir/$sraid.fasta -dbtype nucl`;
		}
	}
	# create 6-frame translation database if needed
	if ( $doProfile == 1 ){
		#if ( ! -e "$tmpdir/$sraid-aa.fasta" ){
		#	printf "[virushunter] \t$sraid: creating 6-frame translation DB for HMMer\n";
		#	`transeq $tmpdir/$sraid.fasta $tmpdir/$sraid-aa.fasta -frame=6`;
		#}
		if ( ! -e "$tmpdir/$sraid-aa-F1.fasta" ){
			printf "[virushunter] \t$sraid: creating 6 single-frame translation DBs for HMMer\n";
			my $pm = new Parallel::ForkManager( 6 );
			for ( my $fi=1; $fi<=6; $fi++ ){
				# begin fork
				my $pid   = $pm->start and next;
				# do the work
				my $frame = $fi;
				$frame = ($frame - 3) * -1  if ( $frame > 3 );
				`$transeqbin $tmpdir/$sraid.fasta $tmpdir/$sraid-aa-F$fi.fasta -frame=$frame -sformat1=pearson`;
				# end fork
				$pm->finish;
			}
			# wait for all forks
			$pm->wait_all_children;
		}
	}	
	# preprocess runtime for this run
	my $etime = time();
	$rtimeP   = ($etime - $stime);
} # end prepare

# run blastn against user-defined set of known viruses to filter them out
sub filter {
	printf "[virushunter] \t$sraid: running filtering blastn using %d threads\n", $threads;
	my $cmd  = "blastn -db $tmpdir/$sraid.fasta -query $filterFile -evalue $filterE -num_threads $threads ";
	   $cmd .= "-outfmt '6 sseqid qseqid pident evalue' ";
	open( BLAST, "$cmd |" );
	while ( my $line = <BLAST> ){
		chomp($line);  next if $line eq "";
		my ( $hid,$qid,$ident,$eval ) = split( /\t/,$line );
		if ( $ident >= $filterIdent ){
			$filterHits{$hid} = ();
			$filterHits{$hid}{'query'}       = $qid;
			$filterHits{$hid}{'seqidentity'} = $ident;
			$filterHits{$hid}{'evalue'}      = $eval;
		}		
	}
	printf "[virushunter] \t$sraid: found %d reads that match to a filtering query\n", scalar keys %filterHits;
} # end filter

# report total runtime and size of SRA dataset
sub reportResources{
	# get runtime of data part and SRA run size from log file
	open( DATA, "<$logtmpdir/$sraid-data.txt" );
	while ( my $line = <DATA> ){
		chomp($line);  next if $line eq "";
		if ( $line =~ /Gb: (.*)/ ){
			$sraSize = $1; # size of (.sra compressed) SRA run
		}
		if ( $line =~ /runtime: (\d+) h (\d+) min (\d+) sec/ ){
			$rtimeD = $1 * 60 * 60 + $2 * 60 + $3; # runtime for data processing
		}
	}
	close(DATA);
	# report runtimes of different parts
	printf "\n[virushunter] \t$sraid: runtime download   part : %s\n", runtime2hms( $rtimeD );
	printf   "[virushunter] \t$sraid: runtime preprocess part : %s\n", runtime2hms( $rtimeP );
	printf   "[virushunter] \t$sraid: runtime search     part : %s\n", runtime2hms( $rtimeS );
	printf   "[virushunter] \t$sraid: runtime filter     part : %s\n", runtime2hms( $rtimeF );
	if ( $mapHg38 ){
	printf   "[virushunter] \t$sraid: runtime mapping    part : %s\n", runtime2hms( $rtimeM );
	}
	printf   "[virushunter] \t$sraid: total runtime           : %s for %.2f Gb\n", runtime2hms( $rtimeD+$rtimeP+$rtimeS+$rtimeF+$rtimeM ), $sraSize;		
} # end reportResources

# runtime for h-min-sec format
sub runtime2hms{
	my $rtime = shift;
	my $ho    = floor( $rtime / 3600 );  $rtime -= $ho * 3600;
	my $mi    = floor( $rtime /   60 );  $rtime -= $mi *   60;
	my $hms   = sprintf "%2d h %2d min %2d sec", $ho, $mi, $rtime;
	return( $hms );
} # end runtime2hms

# run and parse tblastn search
sub searchByBlast {
	# verify that SRA dataset was successfully downloaded and pre-processed
	exit() if ! -e "$tmpdir/$sraid.fasta";
	# init new hits
	my %hits  = ();
	# run Blast
	printf "[virushunter] \t$sraid: running tblastn using %d threads\n", $threads;
	my ($qidN, $qidO, $outp, $alns) = ("", "x", "", "");
	my ($hid, $ident, $alen,$score, $eval, $qseq, $hseq);
	my %hitsPerRead = ();  my %hitseqs = ();  my %hitsQuery = ();
	my $evalMin = 1000000;
	my $hitsIdent = 0;
	my $Ecut = $blastE > 10 ? $blastE : 10;
	my $cmd  = "tblastn -db $tmpdir/$sraid.fasta -query $queryF -evalue $Ecut -num_threads $threads ";
	   $cmd .= "-outfmt '6 qseqid sseqid pident length evalue bitscore qseq sseq' ";
	open( BLAST, "$cmd |" );
	while ( my $line = <BLAST> ){
		chomp($line);  next if $line eq "";
		my ($qidN,$hid,$ident,$alen,$eval,$score,$qseq,$hseq) = split( /\t/,$line );
		# do not consider hits that were found during (optional) filtering step
		next if exists $filterHits{ $hid };
		# save hit sequence for motif search
		my $hseq2 = $hseq;  $hseq2 =~ s/-//g;
		if ( ! exists $hitseqs{$hid} ){
			if ( $alen >= $readHitLenTH ){
				$hitseqs{$hid} = "XXX".$hseq2."XXX";	# add terminal XXX to account for partial motifs
			}
		}elsif( length($hseq) > length($hitseqs{$hid}) ){
			if ( $alen >= $readHitLenTH ){
				$hitseqs{$hid} = "XXX".$hseq2."XXX";	# add terminal XXX to account for partial motifs
			}
		}
		# save globally best E-value
		if ( $eval < $evalMin ){  $evalMin = $eval; }
		# save hit if E-value below threshold
		if ( $eval < $blastE  or  $doRefSeqBlast == 1 ){
			if ( (exists $hits{$hid} and $score > $hits{$hid}{'score'}) or !exists $hits{$hid} ){
				$hits{$hid}{'query'} = $qidN;
				$hits{$hid}{'score'} = $score;
				$hits{$hid}{'eval'}  = $eval;
			}
		}
		# save hits per read and hits with good identity if E-value below threshold and hit length above threshold
		if ( ( $eval < $blastE  or  $doRefSeqBlast == 1 )  and  $alen >= $readHitLenTH ){
			if ( ! exists $hitsPerRead{$hid} ){       $hitsPerRead{$hid} = 0;  }
			if ( ! exists $hitsQuery{$hid} ){         $hitsQuery{$hid}   = (); }
			if ( ! exists $hitsQuery{$hid}{$qidN} ){  $hitsPerRead{$hid}++;    }
			$hitsQuery{$hid}{$qidN} = 1;
			if ( $ident >= $identTH ){  $hitsIdent++; }
		}
		# build custom blast output
		if ( $qidN ne $qidO ){
			if ( $qidO ne "x" ){
				$outp .= sprintf "%-20s\t\t%8s\t%8s\t%8s\t%8s\n\n", "-----","-------","--------","-----","-------";	
				$outp .= $alns."\n";
			}
			$alns  = "";
			$outp .= "============\n";
			$outp .= sprintf "Query= %s\n", $qidN;
			$outp .= "============\n\n";
			$outp .= sprintf "%-20s\t\t%8s\t%8s\t%8s\t%8s\n", "hitID","len[aa]","ident[%]","score","E-value";
			$outp .= sprintf "%-20s\t\t%8s\t%8s\t%8s\t%8s\n", "-----","-------","--------","-----","-------";	
		}
		$outp .= sprintf "%-20s\t\t%8s\t%8s\t%8s\t%8s\n", $hid,$alen,$ident,$score,$eval;
		my $qidNprint = $qidN;
		if ( length($qidN) > 35 ){  $qidNprint = substr( $qidNprint, 0, 35 ); }
		$alns .= sprintf "%-35s:\t%s\n",   $qidNprint, $qseq;
		my @qre = split(//,$qseq);  my @hre = split(//,$hseq);
		$alns .= sprintf "%-35s \t",  "";
		for ( my $i=0; $i<=$#qre; $i++){
			#if ( $qre[$i] eq $hre[$i] ){  $alns .= "|";
			if ( $qre[$i] eq $hre[$i] ){  $alns .= sprintf "%s", $qre[$i];
			}else{                        $alns .= " "; }
		}
		$alns .= "\n";
		$alns .= sprintf "%-35s:\t%s\n\n", $hid,  $hseq;
		# next hit
		$qidO = $qidN;
	}
	$outp .= sprintf "%-20s\t\t%8s\t%8s\t%8s\t%8s\n\n", "-----","-------","--------","-----","-------";	
	$outp .= $alns;
	close(BLAST);

	# optionally compare hit sequences against motifs
	my $motifhits = scanMotifs( $sraid, \%hitseqs );

	# determine hits and save blast output for them
	my $numhits = scalar keys %hits;
	my $maxhitsPerRead = 0;
	$maxhitsPerRead = (sort {$b <=> $a} values %hitsPerRead)[0] if scalar keys %hitsPerRead > 0;
	#if ( ( $numhits >= $hitsTotalTH and $maxhitsPerRead >= $hitsPerReadTH and $motifhits >= $hitsMotifsTH ) or
	#     ( $evalMin < $blastEsure ) ){
	my $iamhit = "false";
	if ( ( $numhits >= $hitsTotalTH and $motifhits >= $hitsMotifsTH ) or
	     ( $maxhitsPerRead >= $hitsPerReadTH ) or
	     ( $evalMin < $blastEsure ) or 
             ( $hitsIdent >= $identTHn ) ){
		$iamhit = "true";
		if ( ! -d "$tmpdir/$projID/$sraid" ){
			mkdir( "$tmpdir/$projID/$sraid" );
		}
		# save full blast output
		open( BLASTOUT, ">$tmpdir/$projID/$sraid/tblastn.out" );
		print BLASTOUT $outp;
		close(BLASTOUT);
		# save hits separately
		my $readsavedN = 0;
		open( BLASTHITS, ">$tmpdir/$projID/$sraid/tblastn-hits.tsv" );
		foreach ( sort {$hits{$a}{'eval'} <=> $hits{$b}{'eval'}} keys %hits ){
			printf BLASTHITS "%s\t%s\t%s\t%s\n", $_,$hits{$_}{'query'},
				                             $hits{$_}{'score'},$hits{$_}{'eval'};
			$readsavedN++;
			last if $readsavedN >= $maxHitN;
		}
		close(BLASTHITS);
		# save info on BioSample and potential other run IDs
		#`echo $sampleID > $resdir/$projID/hits/$sraid/sampleID.txt`;
	}else{
		`rm -rf $tmpdir/$projID/$sraid`;
	}

	# search runtime for this run
	my $etime = time();
	$rtimeS   = ($etime - $stime - $rtimeP);
	
	# report info for this run
	printf "[virushunter] \t$sraid: tblastn hits:              %d\n",   $numhits;
	printf "[virushunter] \t$sraid: max queries hit per read:  %d\n",   $maxhitsPerRead;
	printf "[virushunter] \t$sraid: best E-value:              %.1e\n", $evalMin;
	printf "[virushunter] \t$sraid: hits with >%.1f%% identity: %d\n",  $identTH, $hitsIdent;
	printf "[virushunter] \t$sraid: FIMO motif hits:           %d\n",   $motifhits;
	printf "[virushunter] \t$sraid: detected as hit:           %s\n",   $iamhit;	
} # end searchByBlast

# run and parse HMMer search
sub searchByHMMer {
	# verify that SRA dataset was successfully downloaded
	exit() if ! -e "$tmpdir/$sraid.fasta";
	# init new hits and runtime variables
	my %hits  = ();
	# run HMMer
	#printf "[virushunter] \t$sraid: running hmmsearch using %d threads\n", $threads;
	#my $Ecut = $hmmerE > 10 ? $hmmerE : 10;
	#my $cmd  = "hmmsearch -E $hmmerE --cpu $threads --tblout $tmpdir/$sraid-hmmsearch.tsv $queryF $tmpdir/$sraid-aa.fasta";
	#   $cmd .= " > $tmpdir/$sraid-hmmsearch.txt";
	#`$cmd`;
	#my $threads_frame =  2;
	#   $threads_frame =  4  if ( $threads == 24 );
	my $threads_frame = int( $threads / 6 );
	printf "[virushunter] \t$sraid: running  6 single-frame HMMer searches using %d threads each\n", $threads_frame;
	my $Ecut = $hmmerE > 10 ? $hmmerE : 10;
	my $targetNum  = `grep -c '>' $tmpdir/$sraid-aa-F1.fasta`;  chomp( $targetNum );
	   $targetNum *= 6; # set database size to 6 times the number of reads (because we search each reading frame separately) 
	my $pm = new Parallel::ForkManager( 6 );
	for ( my $fi=1; $fi<=6; $fi++ ){
		# begin fork
		my $pid = $pm->start and next;
		# do the work
		my $cmd  = "hmmsearch -E $hmmerE --cpu $threads_frame --tblout $tmpdir/$sraid-hmmsearch-F$fi.tsv -Z $targetNum $queryF $tmpdir/$sraid-aa-F$fi.fasta";
		   $cmd .= " > $tmpdir/$sraid-hmmsearch-F$fi.txt";
		`$cmd`;
		# end fork
		$pm->finish;
	}
	# wait for all forks
	$pm->wait_all_children;
	# combine results of frame-specific hmmsearches
	`cat $tmpdir/$sraid-hmmsearch-F*.tsv > $tmpdir/$sraid-hmmsearch.tsv`;
	`cat $tmpdir/$sraid-hmmsearch-F*.txt > $tmpdir/$sraid-hmmsearch.txt`;
	# parse HMMer output
	my $evalMin = 1000000;
	open( HMMER, "<$tmpdir/$sraid-hmmsearch.tsv" );
	while ( my $line = <HMMER> ){
		chomp($line);
		next if $line =~ /^#.*/;
		# get hit read id and E-value
		my @v = split( /\s+/, $line );
		$v[0] =~ /(.*)_\d/;
		my $hid  = $1;
		my $eval = $v[4];
		# do not consider hits that were found during (optional) filtering step
		next if exists $filterHits{ $hid };
		# save globally best E-value
		if ( $eval < $evalMin ){  $evalMin = $eval; }
		# save hit if E-value below threshold
		if ( $eval < $hmmerE ){
			$hits{$hid}{'score'}    = $v[5];
			$hits{$hid}{'eval'}     = $eval;
			$hits{$hid}{'qProfile'} = $v[2];
		}
	}
	close(HMMER);
	# determine if this SRA is a hits and save HMMer output for it
	my $numhits = scalar keys %hits;
	my $iamhit = "false";
	if ( ( $numhits >= $hitsTotalHmm ) or
	     ( $evalMin < $hmmerEsure ) ){
		$iamhit = "true";
		if ( ! -d "$tmpdir/$projID/$sraid" ){
			mkdir( "$tmpdir/$projID/$sraid" );
		}
		# save full HMMer output
		`cp $tmpdir/$sraid-hmmsearch.txt $tmpdir/$projID/$sraid/hmmsearch.out`;
		#`cp $tmpdir/$sraid-hmmsearch.txt $resdir/$projID/log/$sraid-hmmsearch.out`;
		# save hits separately
		my $readsavedN = 0;
		open( HITS, ">$tmpdir/$projID/$sraid/hmmsearch-hits.tsv" );
		foreach ( sort {$hits{$a}{'eval'} <=> $hits{$b}{'eval'}} keys %hits ){
			printf HITS "%s\t%s\t%s\t%s\n", $_, $hits{$_}{'qProfile'},
			                                $hits{$_}{'score'}, $hits{$_}{'eval'};
			$readsavedN++;
			last if $readsavedN >= $maxHitN;
		}
		close(HITS);
		# save info on BioSample and potential other run IDs
		#`echo $sampleID > $resdir/$projID/hits/$sraid/sampleID.txt`;
	}else{
		`rm -rf $tmpdir/$projID/$sraid`;
	}
	
	# search runtime for this run
	my $etime = time();
	$rtimeS   = ($etime - $stime - $rtimeP);
	
	# report final info for this run
	printf "[virushunter] \t$sraid: hmmsearch hits:            %d\n",   $numhits;
	printf "[virushunter] \t$sraid: best E-value:              %.1e\n", $evalMin;
	printf "[virushunter] \t$sraid: detected as hit:           %s\n",   $iamhit;
} # end searchByHMMer

# scan blast hit sequences against user-provided motifs using MEME suite
sub scanMotifs {
	return 1 if $motifF eq "";
	my $sraid = $_[0];
	my %seqs  = %{$_[1]};
	# run FIMO from the MEME suite
	my $wdir  = $tmpdir;
	open( FIMO, ">$wdir/$sraid-fimo-seqs.fas" );
	foreach ( keys %seqs ){
		printf FIMO ">%s\n%s\n", $_, $seqs{$_};
	}
	close(FIMO);
	`fimo --o $wdir/$sraid-fimo-tmp $motifF $wdir/$sraid-fimo-seqs.fas > $wdir/$sraid-fimo.log 2>&1`;
	# count number of FIMO hits and return
	my $num = -1;
	open( FIMO, "<$wdir/$sraid-fimo-tmp/fimo.txt");
	while ( my $line = <FIMO> ){
		chomp($line);
		next if $line eq "";
		if ( $num >= 0){
			my @v = split( /\t/,$line );
			$num++ if $v[7] <= $fimoTHq;
		}else{
			$num++;
		}
	}
	close(FIMO);
	# remove temporary files
	if ( $num >= $hitsMotifsTH ){
		mkdir( "$tmpdir/$projID/$sraid" );
		`cp $wdir/$sraid-fimo-tmp/fimo.txt $tmpdir/$projID/$sraid/fimo.txt`;
		`cp $wdir/$sraid-fimo.log          $tmpdir/$projID/$sraid/fimo.log`;
	}
	if ( $debugMode == 0 ){
		`rm -rf $wdir/$sraid-fimo-tmp $wdir/$sraid-fimo-seqs.fas $wdir/$sraid-fimo.log`;
	}
	# return number of found motif occurrences
	return $num;
} # end scanMotifs


# step-wise filtering of hit reads using RefSeq
sub searchRefSeq{
	# only run for SRA experiments that are hits
	return 0 if ! -d "$tmpdir/$projID/$sraid";
	# verify that SRA dataset was successfully downloaded
	return 0 if ! -e "$tmpdir/$sraid.fasta";
	printf "\n[virushunter] \t$sraid: comparing hits against host and viral reference sequences\n";
        # read hits and info from search result file
        my $rawHitsFile = "$tmpdir/$projID/$sraid/tblastn-hits.tsv";
	   $rawHitsFile = "$tmpdir/$projID/$sraid/hmmsearch-hits.tsv" if $doProfile == 1;
	my %readHits = ();
	open( RH, "<$rawHitsFile" ) or die( "Can't open file '$rawHitsFile': $!\n" );
	while ( my $line = <RH> ){
		chomp($line);  next if $line eq "";
		my @v = split( /\t/, $line );
		$readHits{ $v[0] } = ();
		$readHits{ $v[0] }{ 'query' } = $v[1];
		$readHits{ $v[0] }{ 'E' }     = $v[3];
	}
	close(RH);
	# extract reads
	my $rawHitsIDs = "$tmpdir/$sraid-HitReadsAll.ids";
	#my $rawHitsFas0 = "$tmpdir/$sraid-HitReadsAll.fasta";
	my $rawHitsFas = "$tmpdir/$sraid-HitReadsAll.fasta";
	open( RHIDS, ">$rawHitsIDs" ) or die( "Can't open file '$rawHitsIDs': $!\n" );
	foreach ( sort keys %readHits ){  printf RHIDS "%s\n", $_; }
	close( RHIDS);
	my $cmd1 = "grep -F --no-group-separator -A 1 -f $rawHitsIDs $tmpdir/$sraid.fasta > $rawHitsFas";
	`$cmd1`;
	# trim adapters, other contaminant sequences and low-quality bases
	#my $rawHitsFas  = "$tmpdir/$sraid-HitReadsAll-trimmed.fasta";
	#trimReads( $rawHitsFas0, $rawHitsFas );
	#$rawHitsFas = $rawHitsFas0 if ( -s $rawHitsFas == 0);	# return to original untrimmed file if nothing remained after cutadapt
	# assemble reads to contigs
	printf "[virushunter] \t$sraid: assembling initial hits with $assembler\n";
	my $assLog      = $rawHitsFas.".cap.log";
	my $assContigs  = $rawHitsFas.".cap.contigs";
	my $assSinglets = $rawHitsFas.".cap.singlets";
	my $cmd2 = "$assembler $rawHitsFas -h $cap3overhang -o $cap3overlap > $assLog";		
	`$cmd2`;
	# read contigs and reformat
	my %contigs  = ();
	my $contigID = "";
	open( CNTG, "<$assContigs" ) or die ( "Can't open file '$assContigs': $!\n" );
	while ( my $line = <CNTG> ){
		chomp($line);  next if $line eq "";
		if ( $line =~ />(.*)/ ){
			$contigID = $sraid.".".$1.".1";
			$contigs{ $contigID } = {};
			$contigs{ $contigID }->{ 'seq' }       = "";
			$contigs{ $contigID }->{ 'reads' }     = [];
			$contigs{ $contigID }->{ 'E_refseq' }  = 1e6;		# init here, to be able to switch filters   
			$contigs{ $contigID }->{ 'gi_refseq' } = $contigID;	# init here, to be able to switch filters   
		}else{
			$contigs{ $contigID }->{ 'seq' }      .= $line;
		}
	}
	close(CNTG);
	# get read-contig mapping
	open( CAP, "<$assLog" ) or die( "Can't open file '$assLog': $!\n" );
	$contigID = "";
	while ( my $line = <CAP> ){
		chomp($line);  next if $line eq "";
		last if $line eq "DETAILED DISPLAY OF CONTIGS";
		if ( $line =~ /\*+ Contig (\d+) \*+/ ){
			$contigID = $sraid.".Contig".$1.".1";
		}elsif ( $contigID ne "" ){
			$line =~ s/^\s+//; $line =~ s/ .*//;  $line =~ s/\-$//;  $line =~ s/\+$//;
			push( @{ $contigs{ $contigID }->{ 'reads' } }, $line );
		}		
	}
	close(CAP);
	# read singlets and combine with contigs
	open( SNGLT, "<$assSinglets" ) or die ( "Can't open file '$assSinglets': $!\n" );
	while ( my $line = <SNGLT> ){
		chomp($line);  next if $line eq "";
		if ( $line =~ />(.*)/ ){
			$contigID = $1;  $contigID =~ s/ .*//;
			$contigs{ $contigID }->{ 'seq' }       = "";
			$contigs{ $contigID }->{ 'reads' }     = [];
			$contigs{ $contigID }->{ 'E_refseq' }  = 1e6;		# init here, to be able to switch filters   
			$contigs{ $contigID }->{ 'gi_refseq' } = $contigID;	# init here, to be able to switch filters   
			push( @{ $contigs{ $contigID }->{ 'reads' } }, $contigID );
		}else{
			$contigs{ $contigID }->{ 'seq' }      .= $line;
		}
	}
	close(SNGLT);
	# write reformatted contigs and singlets to file
	my $contigIDsRaw = "$tmpdir/$projID/$sraid/contigs.singlets-unfiltered.txt";
	open( CONTSING, ">$contigIDsRaw" ) or die ( "Can't open file '$contigIDsRaw': $!\n" );
	foreach ( sort keys %contigs ){
		printf CONTSING "%s\n", $_;
	}
	close(CONTSING);
	# run filter 1 - blastx against custom contaminant DB
	my $contigIDsFilter1  = "$tmpdir/$projID/$sraid/contigs.singlets-filter1.txt";
	my $contigSeqsFilter1 = "$tmpdir/$sraid-HitReadsAll-trimmed.fasta.cap.contigs.singlets-filter1.fas";
	filter1( $contigIDsRaw,     $contigIDsFilter1, $contigSeqsFilter1, \%contigs );
	# run filter 2 - tblastx against viral-genomic
	my $contigIDsFilter2  = "$tmpdir/$projID/$sraid/contigs.singlets-filter2.txt";
	my $contigSeqsFilter2 = "$tmpdir/$sraid-HitReadsAll-trimmed.fasta.cap.contigs.singlets-filter2.fas";
	#filter2( $contigIDsFilter1, $contigIDsFilter2, $contigSeqsFilter2, \%contigs );
	filter3( $contigIDsFilter1, $contigIDsFilter2, $contigSeqsFilter2, \%contigs );
	# run filter 3 - blastx against non-viral refseq_protein
	my $contigIDsFilter3  = "$tmpdir/$projID/$sraid/contigs.singlets-filter3.txt";
	my $contigSeqsFilter3 = "$tmpdir/$sraid-HitReadsAll-trimmed.fasta.cap.contigs.singlets-filter3.fas";
	#filter3( $contigIDsFilter2, $contigIDsFilter3, $contigSeqsFilter3, \%contigs );
	filter2( $contigIDsFilter2, $contigIDsFilter3, $contigSeqsFilter3, \%contigs );
	# extract remaining hits to file
	my $contigSeqsFinal   = "$tmpdir/$projID/$sraid/contigs.singlets.fas";
	my $contigReadsFinal  = "$tmpdir/$projID/$sraid/contigs.singlets.reads.tsv";
	my @idsFinal = ();
	open( SOUT, ">$contigSeqsFinal" )  or die ( "Can't open file '$contigSeqsFinal': $!\n" );	
	open( ROUT, ">$contigReadsFinal" ) or die ( "Can't open file '$contigReadsFinal': $!\n" );	
	open(  IIN, "<$contigIDsFilter3" ) or die ( "Can't open file '$contigIDsFilter3': $!\n" );
	while ( my $line = <IIN> ){
		chomp($line);  next if $line eq "";
		push( @idsFinal, $line );
		printf SOUT ">%s\n%s\n", $line, $contigs{ $line }->{ 'seq' };
		foreach ( @{ $contigs{ $line }->{ 'reads' } } ){
			printf ROUT "%s\t%s\n", $line, $_;
		}
	}	
	close( IIN );
	close(ROUT );
	close(SOUT );
	# remove this SRA from hit directory if no hits remain after filtering
	if ( scalar @idsFinal == 0 ){
		# filter runtime for this run
		my $etime = time();
		$rtimeF   = ($etime - $stime - $rtimeP - $rtimeS);
		# end here
		`rm -rf $tmpdir/$projID/$sraid`;
		return 0;
	}
	# summarize final hits
	my $res = {};
	my $contig2query = {};
	foreach my $sid ( @idsFinal ){
		# sort by viral refseq subject that was hit
		my $vid = $contigs{ $sid }->{ 'gi_refseq' };
		if ( ! exists $res->{ $vid } ){
			$res->{ $vid } = {};
			$res->{ $vid }->{ 'init_best_E' }       = 1e6;
			$res->{ $vid }->{ 'init_best_query' }   = "";
			$res->{ $vid }->{ 'init_num_hits' }     = 0;
			$res->{ $vid }->{ 'refseq_contigs' }    = 0;
			$res->{ $vid }->{ 'refseq_best_E' }     = 1e6;
			$res->{ $vid }->{ 'refseq_best_Ident' } = 0;
			$res->{ $vid }->{ 'refseq_subject' }    = sprintf "gi:%d|%s", $vid, $contigs{ $sid }->{ 'id_refseq' };
		}
		# best E, query, and number of hit reads of initial search
		$contig2query->{ $sid } = {};
		$contig2query->{ $sid }->{ 'E' }     = 1e6;
		$contig2query->{ $sid }->{ 'query' } = "";
		foreach my $rid ( @{ $contigs{ $sid }->{ 'reads' } } ){
			if ( $readHits{ $rid }{ 'E' } < $res->{ $vid }->{ 'init_best_E' } ){
				$res->{ $vid }->{ 'init_best_E' }     = $readHits{ $rid }{ 'E' };
				$res->{ $vid }->{ 'init_best_query' } = $readHits{ $rid }{ 'query' };
			}
			$res->{ $vid }->{ 'init_num_hits' }++;
			if ( $readHits{ $rid }{ 'E' } < $contig2query->{ $sid }->{ 'E' } ){
				$contig2query->{ $sid }->{ 'E' }     = $readHits{ $rid }{ 'E' };
				$contig2query->{ $sid }->{ 'query' } = $readHits{ $rid }{ 'query' };
			}
		}
		# best E, subject, number of reads and contigs, and subject name of RefSeq search
		if ( $contigs{ $sid }->{ 'E_refseq' } < $res->{ $vid }->{ 'refseq_best_E' } ){
			$res->{ $vid }->{ 'refseq_best_E' }     = $contigs{ $sid }->{ 'E_refseq' };
			$res->{ $vid }->{ 'refseq_best_Ident' } = $contigs{ $sid }->{ 'ident_refseq' };
			$res->{ $vid }->{ 'refseq_best_lens' }  = $contigs{ $sid }->{ 'lens_refseq' };
		}
		$res->{ $vid }->{ 'refseq_contigs' }++;
	}
	# annotate final hits with taxonomy info
	printf "[virushunter] \t$sraid: adding taxonomy annotation for final viral reference sequence hits\n";	
	my $giFile = "$tmpdir/$sraid-FinalHits-gis.txt";
	open( GI, ">$giFile" ) or die( "Can't open file '$giFile': $!\n" );
	foreach my $gi ( keys %{$res} ){
		printf GI "|%d|\n", $gi;
	}
	close(GI);
	my $cmdTax = "$gi2taxExec $giFile";
	open( TAX, "$cmdTax |" );
	while ( my $line = <TAX> ){
		chomp($line);  next if $line eq "";
		my @v = split( /\t/, $line );
		$res->{ $v[0] }->{ 'refseq_taxonomy' } = $v[1];
	}
	close(TAX);
	# write contig to initial_query mapping to file
	my $contig2queryFile = "$tmpdir/$projID/$sraid/contigs.singlets.query.tsv";
	open( C2Q, ">$contig2queryFile" ) or die( "Can't open file '$contig2queryFile': $!\n" );
	foreach my $sid ( sort { $contig2query->{ $a }->{ 'query' } cmp $contig2query->{ $b }->{ 'query' } } keys %{ $contig2query } ){
		printf C2Q "%s\t%s\t%.2e\n", $sid, $contig2query->{ $sid }->{ 'query' }, $contig2query->{ $sid }->{ 'E' } ;
	}
	close(C2Q);	
	# write result summary to file
	my $finalResFile = "$tmpdir/$projID/$sraid/final.hits.tsv";
	open( RES, ">$finalResFile" ) or die( "Can't open file '$finalResFile': $!\n" );
	foreach my $vid ( sort { $res->{$a}->{'init_best_E'} <=> $res->{$b}->{'init_best_E'} } keys %{$res} ){
		printf RES "%d\t%.2e\t%s\t%.2e\t%.1f\t%s\t%d\t%s\t", $res->{$vid}->{'init_num_hits'},  $res->{$vid}->{'init_best_E'},       $res->{$vid}->{'init_best_query'},
							  	     $res->{$vid}->{'refseq_best_E'},  $res->{$vid}->{'refseq_best_Ident'}, $res->{$vid}->{'refseq_best_lens'},
								     $res->{$vid}->{'refseq_contigs'}, $res->{$vid}->{'refseq_subject'};
		printf RES "%s",                                     $res->{$vid}->{'refseq_taxonomy'}  if exists $res->{$vid}->{'refseq_taxonomy'};
		printf RES "\n";
	}
	close(RES);
	# filter runtime for this run
	my $etime = time();
	$rtimeF   = ($etime - $stime - $rtimeP - $rtimeS);
}

# filter - blastx against custom contaminant DB
sub filter1{
	my $idsin  = shift;
	my $idsout = shift;
	my $fasout = shift;
	my $cntgs  = shift;
	#printf "[virushunter] \t$sraid: filtering - tblastx against custom contaminant DB\n";
	printf "[virushunter] \t$sraid: filtering -  blastn against custom contaminant DB\n";
	# read query sequence IDs from file
	open( IDSIN, "<$idsin" ) or die ( "Can't open file '$idsin': $!\n" );
	my %ids = ();
	while ( my $line = <IDSIN> ){
		chomp($line);  next if $line eq "";  $ids{ $line } = 1;
	}
	close(IDSIN);
	# write query sequences to file
	open( FASOUT, ">$fasout" ) or die ( "Can't open file '$fasout': $!\n" );
	foreach ( sort keys %ids ){
		printf FASOUT ">%s\n%s\n", $_, $$cntgs{ $_ }->{ 'seq' };
	}
	close(FASOUT);
	# blast
	#my $blastfile = $fasout."-tblastx.tsv";
	my $blastfile = $fasout."-blastn.tsv";
	#my $cmd = "tblastx -db $filterDB -query $fasout -outfmt '6 qseqid sseqid pident evalue' -evalue $DBblastEh -num_threads $threads > $blastfile";
	my $cmd = "blastn -db $filterDB -query $fasout -outfmt '6 qseqid sseqid pident evalue' -evalue $DBblastEh -num_threads $threads > $blastfile";
	`$cmd`;
	# get hits
	open( BLAST, "<$blastfile" ) or die( "Can't open file '$blastfile': $!\n" );
	while ( my $line = <BLAST> ){
		chomp($line);  next if $line eq "";
		my @v = split( /\t/, $line );
		delete $ids{ $v[0] };
	}
	close(BLAST);
	# write positive hit sequences to file
	open( IDSOUT, ">$idsout" ) or die( "Can't open file '$idsout': $!\n" );
	foreach ( sort keys %ids ){
		printf IDSOUT "%s\n", $_;
	}
	close(IDSOUT);
	# report remaining number of sequences
	printf "[virushunter] \t$sraid: filtering -  %d hits remaining\n", scalar keys %ids;
} # end filter1

# filter - tblastx against viral_genomic
sub filter2{
	my $idsin  = shift;
	my $idsout = shift;
	my $fasout = shift;
	my $cntgs  = shift;
	printf "[virushunter] \t$sraid: filtering - tblastx against viral_genomic\n";
	# read query sequence IDs from file
	open( IDSIN, "<$idsin" ) or die ( "Can't open file '$idsin': $!\n" );
	my %ids = ();
	while ( my $line = <IDSIN> ){
		chomp($line);  next if $line eq "";  $ids{ $line } = 1;
	}
	close(IDSIN);
	# write query sequences to file
	open( FASOUT, ">$fasout" ) or die ( "Can't open file '$fasout': $!\n" );
	foreach ( sort keys %ids ){
		printf FASOUT ">%s\n%s\n", $_, $$cntgs{ $_ }->{ 'seq' };
	}
	close(FASOUT);
	# blast
	my $blastfile = $fasout."-tblastx.tsv";
	my $cmd = "tblastx -db $viralDB -query $fasout -outfmt '6 qseqid sgi stitle pident evalue length qlen sseqid' -max_hsps 1 -max_target_seqs 1 -evalue $DBblastEv -num_threads $threads > $blastfile";
	`$cmd`;
	# get hits
	my %idsO = ();
	open( BLAST, "<$blastfile" ) or die( "Can't open file '$blastfile': $!\n" );
	while ( my $line = <BLAST> ){
		chomp($line);  next if $line eq "";
		my @v = split( /\t/, $line );
		$idsO{   $v[0] } = 1;
		$$cntgs{ $v[0] }->{ 'E_refseq' }     = $v[4];
		$$cntgs{ $v[0] }->{ 'ident_refseq' } = $v[3];
		$$cntgs{ $v[0] }->{ 'id_refseq' }    = $v[2];
		$$cntgs{ $v[0] }->{ 'id_refseq' }    = $v[7]  if $v[2] eq "N/A";		
		$$cntgs{ $v[0] }->{ 'gi_refseq' }    = $v[1];
		$$cntgs{ $v[0] }->{ 'lens_refseq' }  = sprintf "%d / %d", $v[5]*3, $v[6];
	}
	close(BLAST);
	# save hits
	#print "cp $blastfile $tmpdir/$projID/$sraid/tblastx-viral_genomic-hits.tsv";
	`cp $blastfile $tmpdir/$projID/$sraid/tblastx-viral_genomic-hits.tsv`;
	# write positive hit sequences to file
	open( IDSOUT, ">$idsout" ) or die( "Can't open file '$idsout': $!\n" );
	foreach ( sort keys %idsO ){
		printf IDSOUT "%s\n", $_;
	}
	close(IDSOUT);
	# report remaining number of sequences
	printf "[virushunter] \t$sraid: filtering -  %d hits remaining\n", scalar keys %idsO;
} # end filter2

# filter - blastx against non-viral refseq_protein
sub filter3{
	my $idsin  = shift;
	my $idsout = shift;
	my $fasout = shift;
	my $cntgs  = shift;
	printf "[virushunter] \t$sraid: filtering -  blastx against non-viral refseq_protein\n";
	# read query sequence IDs from file
	open( IDSIN, "<$idsin" ) or die ( "Can't open file '$idsin': $!\n" );
	my %ids = ();
	while ( my $line = <IDSIN> ){
		chomp($line);  next if $line eq "";  $ids{ $line } = 1;
	}
	close(IDSIN);
	# only analyze contig with best E value per refseq subject
	my %ids2 = ();
	foreach ( keys %{$cntgs} ){
		my $vid = $$cntgs{ $_ }->{ 'gi_refseq' };
		if ( ! exists $ids2{ $vid } ){
			$ids2{ $vid } = ();
			$ids2{ $vid }{ 'contigID' } = $_;
			$ids2{ $vid }{ 'Evalue' }   = $$cntgs{ $_ }->{ 'E_refseq' };
		}elsif( $$cntgs{ $_ }->{ 'E_refseq' } < $ids2{ $vid }{ 'Evalue' } ){
			$ids2{ $vid }{ 'contigID' } = $_;
			$ids2{ $vid }{ 'Evalue' }   = $$cntgs{ $_ }->{ 'E_refseq' };
		}
	}
	# write query sequences to file
	open( FASOUT, ">$fasout" ) or die ( "Can't open file '$fasout': $!\n" );
	foreach ( sort keys %ids2 ){
		my $cid = $ids2{ $_ }{ 'contigID' };
		printf FASOUT ">%s\n%s\n", $cid, $$cntgs{ $cid }->{ 'seq' };
	}
	close(FASOUT);
	# blast
	my $blastfile = $fasout."-blastx.tsv";
	my $cmd  = "blastx -db $refseqDB -query $fasout -outfmt '6 qseqid sacc stitle pident evalue' -max_hsps 1 -max_target_seqs 1 -evalue $DBblastEh -num_threads $threads ";
	   $cmd .= "-negative_gilist $refseqNegGIs > $blastfile";
	`$cmd`;
	# get hits
	open( BLAST, "<$blastfile" ) or die( "Can't open file '$blastfile': $!\n" );
	while ( my $line = <BLAST> ){
		chomp($line);  next if $line eq "";
		my @v = split( /\t/, $line );
		next if $v[2] =~ /.*virus.*/i;
		my $cid = $v[0];
		my $vid = $$cntgs{ $cid }->{ 'gi_refseq' };
		# delete all contigs that hit the respective refseq subject
		foreach ( keys %{$cntgs} ){
			if ( $$cntgs{ $_ }->{ 'gi_refseq' } eq $vid ){
				delete $ids{ $_ };
			}
		}
	}
	close(BLAST);
	# write positive hit sequences to file
	open( IDSOUT, ">$idsout" ) or die( "Can't open file '$idsout': $!\n" );
	foreach ( sort keys %ids ){
		printf IDSOUT "%s\n", $_;
	}
	close(IDSOUT);
	# report remaining number of sequences
	printf "[virushunter] \t$sraid: filtering -  %d hits remaining\n", scalar keys %ids;
} # end filter3

# trim reads using autoadapt
sub trimReads{
	my $infile  = shift;
	my $outfile = shift;
	printf "[virushunter] \t$sraid: trimming adapter sequences and low-quality bases for initial hits\n";
	# prepare formatted read IDs
	my $readIDfile = "$infile.fqids";
	`grep '>' $infile | sed 's/>/@/' > $readIDfile`;
	# temporary extract fastq reads
	my $fastqFile1 = "$infile";
	   $fastqFile1 =~ s/\.fasta/\.fastq/;
	`fastq-dump -B --split-spot --skip-technical --readids -O $tmpdir $tmpdir/$sraid.sra`;
	`grep -F -f $readIDfile --no-group-separator -A 3 $tmpdir/$sraid.fastq > $fastqFile1`;
	`rm $tmpdir/$sraid.fastq`;
	# run autoadapt
	my $fastqFile2   = "$fastqFile1";
	   $fastqFile2   =~ s/\.fastq/-trimmed\.fastq/;
	my $autoadaptLog = "$infile.autoadapt.log";
	my $wdir = cwd;
	chdir( $tmpdir );
	`autoadapt --threads=$threads --quality-cutoff=$cutQth $fastqFile1 $fastqFile2 > $autoadaptLog`;
	chdir( $wdir );
	# convert fastq to fasta
	`seqtk seq -a $fastqFile2 > $outfile`;
} # end trimReads

# extract top read hit(s) and search against RefSeq proteins
sub searchRefSeq_old{
	# only run for SRA experiments that are hits
	return 0 if ! -d "$tmpdir/$projID/$sraid";	
	# verify that SRA dataset was successfully downloaded
	return 0 if ! -e "$tmpdir/$sraid.fasta";
	printf "\n[virushunter] \t$sraid: comparing hits against host and viral reference sequences\n";	
	my $rFile1 = "$tmpdir/$sraid-HitReadsAll.fasta";
	if ( $doProfile == 1){
		# get ID of all read hit
		my $hitFile    = "$tmpdir/$projID/$sraid/hmmsearch-hits.tsv";
		my $readidFile = "$tmpdir/$projID/$sraid/hmmsearch.ids";
		open( H, "<$hitFile" )    or die( "Can't open file '$hitFile': $!\n" );
		open( R, ">$readidFile" ) or die( "Can't open file '$readidFile': $!\n" );
		my $line = <H>;  chomp($line);
		my @line = split( "\t", $line );
		my $rid  = $line[0];
		print R $rid."\n";
		while ( $line = <H> ){
			chomp($line);  next if $line eq "";	
			@line = split( "\t", $line );
			print R $line[0]."\n";
		}
		close(R);
		close(H);
		# get read sequence to be used in RefSeq search
		my $cmd = "grep -F --no-group-separator -A 1 -f $readidFile $tmpdir/$sraid.fasta > $rFile1";
		`$cmd`;
	}else{
		# get ID of all read hits for each query
		my $hitFile = "$tmpdir/$projID/$sraid/tblastn.out";
		open( H, "<$hitFile" ) or die( "Can't open file '$hitFile': $!\n");
		my %hreads = ();
		my $qid    = "";
		my $ison   = 0;
		while( my $line = <H> ){
			chomp($line);
			if ( $line =~ /Query= (.*)/ ){              $qid  = $1; $ison = 2; }
			if ( $line =~ /^-----.*/ and $ison == 2 ){  $ison =  1; next; }
			if ( $line =~ /^-----.*/ and $ison == 1 ){  $ison =  0; next; }
			if ( $ison == 1 ){
				my @v   = split( /\s+/, $line );
				my $rid0 = $v[0];
				my $rid  = $rid0;  $rid =~ s/ .*//;
				my $eva  = $v[4];
				if ( ! exists $hreads{$rid} ){
					$hreads{$rid} = ();
					$hreads{$rid}{'Q'} = $qid;
					$hreads{$rid}{'E'} = $eva;
					$hreads{$rid}{'R'} = $rid0;
				}elsif( $hreads{$rid}{'E'} > $eva ){
					$hreads{$rid}{'Q'} = $qid;
					$hreads{$rid}{'E'} = $eva;
					$hreads{$rid}{'R'} = $rid0;
				}
			}
		}
		close(H);
		# remember IDs of all hit reads
		my $readidFile = "$tmpdir/$projID/$sraid/tblastn.ids";
		open( RIDS, ">$readidFile" );
		foreach my $rid ( sort { $hreads{$a}{'E'} <=> $hreads{$b}{'E'} } keys %hreads ){
			print RIDS $rid."\n";
		}
		close(RIDS);
		# get read sequences to be used in RefSeq search
		#my %querydone = ();
		#my $rnum = 0;
		#my $readidFileTmp = "$tmpdir/$sraid-HitReadIDs.tmp";
		#open( RIDS, ">$readidFileTmp" );
		#foreach my $rid ( sort { $hreads{$a}{'E'} <=> $hreads{$b}{'E'} } keys %hreads ){
		#	last if $rnum >= $RefSeqBlastMx;
		#	my $q = $hreads{$rid}{'Q'};
		#	if ( ! exists $querydone{$q} ){
		#		printf RIDS "%s \n", $hreads{$rid}{'R'};
		#		$querydone{$q} = 1;
		#		$rnum++;
		#	}
		#}
		#close(RIDS);
		my $cmd = "grep -F --no-group-separator -A 1 -f $readidFile $tmpdir/$sraid.fasta > $rFile1";
		`$cmd`;
	}
	# run blastx against known (host) proteins to be filtered-out
	printf "[virushunter] \t$sraid: blastx against contaminants\n";	
	my $oFile1 = "$tmpdir/$projID/$sraid/blastx-filter.tsv";
	my $cmdB1  = "blastx -db $filterDB -query $rFile1 -outfmt '6 qseqid sseqid evalue pident' -evalue $DBblastEh -num_threads $threads > $oFile1";
	`$cmdB1`;
	my %filter = ();
	open( FILTER, "<$oFile1" ) or die( "Can't open file '$oFile1': $!\n" );
	while ( my $line = <FILTER> ){
		chomp($line);  next if $line eq "";
		my @v  = split( /\t/, $line );
		my $ri = $v[0]; $ri =~ s/ .*//;
		$filter{$ri} = $v[1];
	}
	close(FILTER);
	# intersect read lists
	my $readidFile2       = "$tmpdir/$projID/$sraid/tblastn.ids";
	my $readidFile2filter = "$tmpdir/$projID/$sraid/tblastn-filtered.ids";
	if ( $doProfile == 1 ){
		$readidFile2       = "$tmpdir/$projID/$sraid/hmmsearch.ids";
		$readidFile2filter = "$tmpdir/$projID/$sraid/hmmsearch-filtered.ids";
	}
	open( RIN,  "<$readidFile2" )       or die( "Can't open file '$readidFile2': $!\n" );
	open( ROUT, ">$readidFile2filter" ) or die( "Can't open file '$readidFile2filter': $!\n" );
	while ( my $line = <RIN> ){
		chomp($line);  next if $line eq "";
		$line =~ / .*/;
		print ROUT $line."\n" if ! exists $filter{$line};
	}
	close(ROUT);
	close(RIN);
	# extract reads to be kept
	my $rFile2    = "$tmpdir/$sraid-HitReadsAll-filtered.fasta";
	my $cmdFilter = "grep -F --no-group-separator -A 1 -f $readidFile2filter $rFile1 > $rFile2";
	`$cmdFilter`;
	# run tblastx against known virus reference genomes
	printf "[virushunter] \t$sraid: tblastx against viral reference genomes\n";	
	my $oFile2 = "$tmpdir/$projID/$sraid/tblastx-ViralGenomic.tsv";
	my $cmdB2  = "tblastx -db $viralDB -query $rFile2 -outfmt '6 evalue pident sgi stitle qseqid' -max_hsps 1 -max_target_seqs 1 -evalue $DBblastEv -num_threads $threads > $oFile2";
	`$cmdB2`;
	# parse tblastx results
	my %tbxRes  = ();
	my %tbxRtoS = ();
	open( TBXI, "<$oFile2" ) or die( "Can't open file '$oFile2': $!\n" );
	while ( my $line = <TBXI> ){
		chomp($line);  next if $line eq "";
		my @v  = split( /\t/, $line );
		my $id = $v[2];
		my $ri = $v[4];  $ri =~ s/ .*//;
		if ( ! exists $tbxRes{$id} ){
			$tbxRes{$id} = ();
			$tbxRes{$id}{'S'} = $v[3];
			$tbxRes{$id}{'N'} = 0;
			$tbxRes{$id}{'E'} = 100;
			$tbxRes{$id}{'I'} = 0;
			$tbxRes{$id}{'R'} = "";
		}
		$tbxRes{$id}{'N'}++;
		$tbxRtoS{$ri}     = $id   if $v[0] < $tbxRes{$id}{'E'};
		$tbxRes{$id}{'R'} = $ri   if $v[0] < $tbxRes{$id}{'E'};
		$tbxRes{$id}{'E'} = $v[0] if $v[0] < $tbxRes{$id}{'E'};
		$tbxRes{$id}{'I'} = $v[1] if $v[1] > $tbxRes{$id}{'I'};
	}
	close(TBXI);
	# run blastx against refseq without viral (only best read per subject)
	printf "[virushunter] \t$sraid: blastx against host reference proteins\n";	
	my $readidFileTop = "$tmpdir/$sraid-HitReadsTop.ids";
	my $rFileTop      = "$tmpdir/$sraid-HitReadsTop.fasta";
	my $oFileTop      = "$tmpdir/$sraid-HitReadsTop-blastx.tsv";
	open( RTOP, ">$readidFileTop" ) or die( "Can't open file '$readidFileTop': $!\n" );
	foreach ( keys %tbxRes ){  printf RTOP "%s\n", $tbxRes{$_}{'R'}; }
	close(RTOP);
	my $cmdTop1  = "grep -F --no-group-separator -A 1 -f $readidFileTop $rFile1 > $rFileTop";
	`$cmdTop1`;
	my $cmdTop2  = "blastx -db $refseqDB -query $rFileTop -outfmt '6 evalue pident sacc stitle qseqid' -max_hsps 1 -max_target_seqs 1 -evalue $DBblastEh -num_threads $threads ";
	   $cmdTop2 .= "-negative_gilist $refseqNegGIs > $oFileTop";
	`$cmdTop2`;
	# filter viral hits using results of refseq blast
	open( BTOP, "<$oFileTop" ) or die( "Can't open file '$oFileTop': $!\n" );
	while ( my $line = <BTOP> ){
		chomp($line);  next if $line eq "";
		my @v  = split( /\t/, $line );
		next if $v[3] =~ /.*virus.*/i;
		my $ri = $v[4];  $ri =~ s/ .*//;
		delete $tbxRes{ $tbxRtoS{$ri} };
	}
	close(BTOP);
	# add taxonomy annotation for hits
	printf "[virushunter] \t$sraid: adding taxonomy annotation for final viral reference sequence hits\n";	
	my $giFile = "$tmpdir/$sraid-HitReadsTop-gis.txt";
	open( GI, ">$giFile" ) or die( "Can't open file '$giFile': $!\n" );
	foreach my $gi ( keys %tbxRes ){
		printf GI "|%d|\n", $gi;
	}
	close(GI);
	my $cmdTax = "$gi2taxExec $giFile";
	open( TAX, "$cmdTax |" );
	while ( my $line = <TAX> ){
		chomp($line);  next if $line eq "";
		my @v = split( /\t/, $line );
		$tbxRes{ $v[0] }{'T'} = $v[1];
	}
	close(TAX);
	# write results to file
	my $oFileFinal      = "$tmpdir/$projID/$sraid/tblastx-ViralGenomic-summary.tsv";
	my $readidFileFinal = "$tmpdir/$projID/$sraid/tblastx-ViralGenomic-reads.ids";
	open( TBXR, ">$readidFileFinal" ) or die( "Can't open file '$readidFileFinal': $!\n" );
	open( TBXO, ">$oFileFinal" )      or die( "Can't open file '$oFileFinal': $!\n" );
	print TBXO "best_E\tbest_Ident\tnum_reads\tsubject\ttaxonomy\n";
	foreach my $id ( sort { $tbxRes{$b}{'I'} <=> $tbxRes{$a}{'I'} } keys %tbxRes ){
		my $cmdRR = "grep -P '\\t$id\\t' $oFile2 | awk 'BEGIN {FS=\"\\t\"}; {print \$5};' ";
		open( CMDRR, "$cmdRR |" ) or die ( "Can't run command '$cmdRR': $!\n" );
		my $hitreadnum = 0;
		while ( my $line = <CMDRR> ){
			chomp($line);  next if $line eq "";
			$hitreadnum++;
			print TBXR "$line\n";
		}
		close(CMDRR);
		#print "[virushunter] \tcmdRR: $cmdRR\n";
		if( $hitreadnum > $tbxRes{$id}{'N'} ){  $tbxRes{$id}{'N'} = $hitreadnum; }
		printf TBXO "%.2e\t%.1f\t%d\tgi:%d|%s\t%s\n", $tbxRes{$id}{'E'}, $tbxRes{$id}{'I'}, $tbxRes{$id}{'N'}, $id, $tbxRes{$id}{'S'}, $tbxRes{$id}{'T'};
	}
	close(TBXO);
	close(TBXR);
	# filter runtime for this run
	my $etime = time();
	$rtimeF   = ($etime - $stime - $rtimeP - $rtimeS);

} # end searchRefSeq_old

# map all reads against reference transcriptome
sub mapRefGenes{
	# only run for SRA experiments that are hits
	return 0 if ! -d "$tmpdir/$projID/$sraid";
	# verify that SRA dataset was successfully downloaded
	return 0 if ! -e "$tmpdir/$sraid.fasta";
	printf "\n[virushunter] \t$sraid: mapping reads against reference transcripts\n";
	# quantify via salmon
	my $cmd  = "salmon quant -l A -i $transcriptsHg38 -r $tmpdir/$sraid.fasta -p $threads";
	   $cmd .= " -g $genemappingHg38 -o $tmpdir/$sraid-hg38-salmon > $tmpdir/$sraid-hg38-salmon.log 2>&1";
	`$cmd`;
	# save result file
	#`cp $tmpdir/$sraid-hg38-salmon/quant.genes.sf $tmpdir/$projID/$sraid/hg38.counts`;
	`awk '{print \$1,\$4,\$5}' $tmpdir/$sraid-hg38-salmon/quant.genes.sf > $tmpdir/$projID/$sraid/hg38.counts`;
	# mapping runtime for this run
	my $etime = time();
	$rtimeM   = ($etime - $stime - $rtimeP - $rtimeS - $rtimeF);
} # end mapRefGenes

# function to compress some result files to save disk space
sub compress{
	# Blast/HMMer output
	if ( -e "$tmpdir/$projID/$sraid/tblastn.out" ){                    `gzip $tmpdir/$projID/$sraid/tblastn.out;   gzip $tmpdir/$projID/$sraid/tblastn-hits.tsv`; }
	if ( -e "$tmpdir/$projID/$sraid/hmmsearch.out" ){                  `gzip $tmpdir/$projID/$sraid/hmmsearch.out; gzip $tmpdir/$projID/$sraid/hmmsearch-hits.tsv`; }
	if ( -e "$tmpdir/$projID/$sraid/tblastx-viral_genomic-hits.tsv" ){ `gzip $tmpdir/$projID/$sraid/tblastx-viral_genomic-hits.tsv`; }
	# other result files
	if ( -e "$tmpdir/$projID/$sraid/contigs.singlets.fas" ){        `gzip $tmpdir/$projID/$sraid/contigs.singlets.fas`; }
	if ( -e "$tmpdir/$projID/$sraid/contigs.singlets.reads.tsv" ){  `gzip $tmpdir/$projID/$sraid/contigs.singlets.reads.tsv`; }
	if ( -e "$tmpdir/$projID/$sraid/contigs.singlets.query.tsv" ){  `gzip $tmpdir/$projID/$sraid/contigs.singlets.query.tsv`; }
	if ( -e "$tmpdir/$projID/$sraid/hg38.counts" ){                 `gzip $tmpdir/$projID/$sraid/hg38.counts`; }
} # end compress

