#!/usr/bin/perl

##########################################################
# assemble viral sequences from NGS data
# workflow optimzed for usage on a computing server
#
# author: Chris Lauber
#  email: chris.lauber@twincore.de
##########################################################

# ------------ #
# load modules
# ------------ #
use warnings;
use strict;
#use Data::Dumper;


# --------- #
# parameter
# --------- #
my $GenSeedExec   = "genseed-hmm";
my $GS_seed_size  = 30;
my $GS_blast_E    = 0.001;
my $GS_blast_ide  = 85;
my $GS_blast_ws   = 7;
my $GS_contig_len = 50000;
my $GS_max_itr    = 100;
my $GS_aln_th     = 75;
my $GS_use_qual   = "no";
my $GS_starter    = 0;
my $GS_direction  = "b";
my $cap_identity  = 85;
my $cap_overhang  = 99; # in %
my $cap_overlap   = 17; # in nt

my $resdir;
my $logdir;
my $LOGF;


# --------------- #
# usage and input
# --------------- #
if ( $#ARGV < 2 ){
die("
usage: virusgathererTWC.pl <parameters>\n
parameters:
\t<family>
\t<fastq_file>
\t<project_id>
\t<run_id>
\t<seed_file>
\t<tmp_dir>
\t<threads>
\t<assembler>
\t<refseqDB>
\t<viralRefSeqIDs>
\t<debugMode> [0|1]
\n");
}


# -------------- #
# get parameters
# -------------- #
my $family        = shift;
my $fastqF        = shift;
my $projID        = shift;
my $sraid         = shift;
my $hunterContigF = shift;
my $tmpdir        = shift;
my $threads       = shift;
my $assembler     = shift;
my $refseqDB      = shift;
my $refseqNegIDs  = shift;
my $debugMode     = shift;


# -------------- #
# initialization
# -------------- #
my ($stime, $etime, $rtime, $ho, $mi);
$stime = time();
init();


# ------------ #
# prepare data
# ------------ #
preprocess();


# --------------------- #
# run iterative assemly
# --------------------- #
runGenseed();
#runTrinity()  if  $TNT == 1;


# ------------------------------------- #
# compare final contigs to reference DB
# ------------------------------------- #
runBlastViral();


# -------- #
# clean-up
# -------- #
cleanup();


# --- #
# end
# --- #
end();


# --------- #
# functions
# --------- #

# run iterative assembly using GenSeed-HMM
sub runGenseed{
	# define seed for Genseed-HMM
	printf "[virusgatherer] using virushunter contigs as seed\n";
	my $seedFile = "$resdir/seed/virushunter.fasta";
	# iteratively assemble seed contigs/singlets to super-contigs
	my ( $capLog, $capContigs, $capSinglets, $capMerged );
	my $capNumOld = 0;
	my $capNumNew = `grep -c '>' $seedFile`;  chomp( $capNumNew );
	my $capItr    = 0;
	while ( $capNumOld != $capNumNew  and  $capItr < 5 ){
		$capLog      = $seedFile.".cap.log";
		$capContigs  = $seedFile.".cap.contigs";
		$capSinglets = $seedFile.".cap.singlets";
		$capMerged   = $seedFile.".cap.contigs.singlets";
		`cap3 $seedFile -h $cap_overhang -o $cap_overlap > $capLog`;
		`cat $capContigs $capSinglets > $capMerged`;
		$capNumOld = $capNumNew;
		$capNumNew = `grep -c '>' $capMerged`;  chomp( $capNumNew );
		$seedFile  = $capMerged;
		$capItr++;
	}
	`mv $seedFile $resdir/seed/seed.fasta`;
	$seedFile = "$resdir/seed/seed.fasta";
	# run GenSeed-HMM
	printf $LOGF "[virusgatherer] running GenSeed-HMM for iterative assembly\n";
	my $cmd  = "$GenSeedExec"; 
	   $cmd .= sprintf " -seed %s",				$seedFile;
	   $cmd .= sprintf " -db %s",				"$resdir/data/$sraid.fasta";
	   $cmd .= sprintf " -assembler %s",			$assembler;
	   $cmd .= sprintf " -ext_seed_size %d",		$GS_seed_size;
	   $cmd .= sprintf " -output %s",			"$resdir/results-$assembler";
	   $cmd .= sprintf " -blastn_parameters \"%s\"",	"-evalue $GS_blast_E -perc_identity $GS_blast_ide -word_size $GS_blast_ws -dust no -soft_masking false";
	   $cmd .= sprintf " -tblastn_parameters \"%s\"",	"-max_target_seqs 50000";
	   $cmd .= sprintf " -hmmsearch_parameters \"%s\"",	"--cpu $threads";
	   $cmd .= sprintf " -clean no";
	   $cmd .= sprintf " -duplicate_headers no";
	   $cmd .= sprintf " -threads %d",			$threads;
	   $cmd .= sprintf " -max_contig_length %d",		$GS_contig_len;
	   $cmd .= sprintf " -max_number_rounds %d",		$GS_max_itr;
	   $cmd .= sprintf " -mapping yes";
	   $cmd .= sprintf " -align_threshold %d",		$GS_aln_th;
	   $cmd .= sprintf " -use_qual %s",			$GS_use_qual;
	   $cmd .= sprintf " -assembler_parameters \"-o $cap_overlap -h $cap_overhang -p $cap_identity\""  if ( $assembler eq "cap3" );
	   $cmd .= sprintf " -aux_starter yes",     if ( $GS_starter == 1 );
	   $cmd .= sprintf " -exp_direction left",  if ( $GS_direction eq "l" );
	   $cmd .= sprintf " -exp_direction right", if ( $GS_direction eq "r" );
	#printf "[virusgatherer] command: %s\n", $cmd;
	`$cmd`;
} # end runGenseed


# assemble full data set using Trinity
# ! note: currently not used !
#sub runTrinity{
#	# prepare files and parameters
#	my $maxMem  = int( $memPerCPU * $threads / 1024 );
#	my @fqFiles = glob( "$tmpdir/data/*trim*.fastq" );
#	my $cmd  = sprintf "trinity --seqType fq --CPU %d --max_memory %dG --output %s", $threads, $maxMem, "$tmpdir/results-trinity";
#	   $cmd .= sprintf " --left %s --right %s", $fqFiles[0], $fqFiles[1]  if $#fqFiles == 1;
#	   $cmd .= sprintf " --single %s",          $fqFiles[0],              if $#fqFiles == 0;
#	# run Trinity
#	printf "[virusgatherer] running Trinity\n";
#	printf "$virushunterTWC.plcmd\n";
#	`$cmd > $tmpdir/trinity.log 2>&1`;
#	# final contigs
#	$finalContigFile = "$tmpdir/results-trinity/Trinity.fasta";
#} # end runTrinity


# blast final contigs against viral reference database
sub runBlastViral{
	# files
	my $finalContigFile = "$resdir/results-$assembler/final_result_dir/final_positive_contigs.fasta";
	my $blastResFile    = "$finalContigFile-vs-viral-blastx.tsv";
	# run blast
	if ( -e $finalContigFile ){
		printf $LOGF "[virusgatherer] Running blast of final contigs against viral reference DB\n";
		my $blastcmd  = "blastx -query $finalContigFile -db $refseqDB -seqidlist $refseqNegIDs";
		   $blastcmd .= " -num_threads $threads -max_hsps 1 -max_target_seqs 1 ";
		   $blastcmd .= " -outfmt '6 qseqid qlen evalue pident length sacc stitle'";
		`$blastcmd > $blastResFile`;
	}
	else{
		printf $LOGF "[virusgatherer] No final contig file to be used by blast found\n";
	}
} # end runBlastViral


# initialization
sub init{
	# create directories not already present
	if( ! -d "$tmpdir/$family/$projID/logs/virusgatherer" ){
		`mkdir $tmpdir/$family/$projID/logs/virusgatherer`;
	}
	if( ! -d "$tmpdir/$family/$projID/results/$sraid/virusgatherer" ){
		`mkdir $tmpdir/$family/$projID/results/$sraid/virusgatherer`;
	}
	if( ! -d "$tmpdir/$family/$projID/results/$sraid/virusgatherer/data" ){
		`mkdir $tmpdir/$family/$projID/results/$sraid/virusgatherer/data`;
	}
	if( ! -d "$tmpdir/$family/$projID/results/$sraid/virusgatherer/seed" ){
		`mkdir $tmpdir/$family/$projID/results/$sraid/virusgatherer/seed`;
	}

	# set directories
	$resdir = "$tmpdir/$family/$projID/results/$sraid/virusgatherer";
	$logdir = "$tmpdir/$family/$projID/logs/virusgatherer";

	# open log file
	open( $LOGF, ">$logdir/$sraid"."_virusgatherer.log" ) or die( "Can't write to log file: $!\n" );

	# record name of machine and data directory
	my $hostname = `hostname`; chomp( $hostname );
	printf $LOGF "[virusgatherer] \t$sraid: running on %s\n", $hostname;
	printf $LOGF "[virusgatherer] \t$sraid: using data directory %s\n\n", $resdir;

} # end init


# prepare fastq data and seed file
sub preprocess{
	# fastq to fasta data transformation incl. trimming
	printf $LOGF "[virusgather] \t$sraid: transforming to Fasta format\n";
	my $read1id = `zcat $fastqF | head -n 1`;
	my $read2id = `zcat $fastqF | head -n 5 | tail -n 1`;
	chomp($read1id);
	chomp($read2id);
	$read1id =~ /(.*\.\d+)\.\d .*/;  my $r1id = $1;
	$read2id =~ /(.*\.\d+)\.\d .*/;  my $r2id = $1;
	my $fastp_w   = $threads;
	   $fastp_w   = 16  if $threads > 16;
	my $fastp_cmd = "fastp -i $fastqF -w $fastp_w -j $resdir/$sraid.fastp.json -h $resdir/$sraid.fastp.html";
	if ( $r1id eq $r2id ){
		$fastp_cmd .= " --interleaved_in --stdout > $resdir/$sraid.trim.fastq";
	}else{
		$fastp_cmd .= " -o $resdir/$sraid.trim.fastq";
	}
	`$fastp_cmd 2>/dev/null`;
	`seqtk seq -A $resdir/$sraid.trim.fastq > $resdir/data/$sraid.fasta 2>/dev/null`;
	`rm $resdir/$sraid*`  if ( $debugMode == 0 );

	# seed file
	`zcat $hunterContigF > $resdir/seed/virushunter.fasta`;
}


# clean up stuff
sub cleanup{
	# delete temporary files
	if ( $debugMode == 0 ){
		`rm -rf $resdir/seed/virushunter.fasta.*`;
	}

	# log end time and calculate runtime
	$etime = time();
	$rtime = $etime - $stime;
	$ho = int( $rtime / 3600 );  $rtime -= $ho * 3600;
	$mi = int( $rtime /   60 );  $rtime -= $mi *   60;
	printf $LOGF "[virusgatherer] done - runtime: %d h %d min %d sec\n", $ho, $mi, $rtime;
} # end cleanup


# do some last things before ending the scirpt
sub end{
	# put empty contig file if not present (needed for snakemake)
	#`touch $resdir/contigs.singlets.fas.gz` if ( ! -e "$resdir/contigs.singlets.fas.gz" );
	# write file showing the search is done
	`echo completed > $resdir/virusgatherer.done`;
	# close log filehandle
	close( $LOGF );
	# exit with code zero
	exit( 0 );
} # end end

