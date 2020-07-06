#!/usr/bin/perl

##########################################################
# iterative assembly using GenSeed-HMM as core engine
# OR
# using Trinity (for RNA-seq data sets)
#
# author: Chris Lauber
##########################################################

# ------------ #
# load modules
# ------------ #
use warnings;
use strict;
#use Data::Dumper;

use lib "/home/lauber/lauber-2015-virusHunter";
use      virusHGutils;


# --------- #
# parameter
# --------- #
my $resdir0       = "/projects/p_sra";
my $tmpdir        = sprintf "%s/lauber-2015-virusGatherer", virusHGutils::get_workspace_path( "scratch" );
my $tmpdir1       = "/tmp";
my $tmpdir2       = sprintf "%s/virusgatherer",             virusHGutils::get_workspace_path( "lustre" );
my $sizeInflatF   = 7; # determined empirically
my $profiledir    = "/home/lauber/lauber-2015-virusHunter/profiles";
my $querydir      = "/home/lauber/lauber-2015-virusHunter/queries";
my $virusdbdir    = "/home/lauber/lauber-2015-virusHunter/virusDB";
my $threads       = 6;
my $useBlastSeed  = 0;
my $useHitReads   = "";
my $userSeed      = "";
my $memPerCPU     = 1875;
my $GenSeedExec   = "genseed-hmm";
my $GS_assembler  = "newbler";
my $GS_seed_size  = 30;
my $GS_blast_E    = 0.001;
my $GS_blast_ide  = 85;
my $GS_blast_ws   = 7;
my $GS_contig_len = 7500;
my $GS_max_itr    = 50;
my $GS_aln_th     = 75;
my $GS_use_qual   = "no";
my $GS_starter    = 0;
my $GS_direction  = "b";
my $cap_identity  = 85;
my $cap_overhang  = 66; # in %
my $cap_overlap   = 17; # in nt
my $TNT           = 0;	# assemble full dataset using Trinity or SPAdes instead of Genseed)
my $useLustre     = 0;
my $userOutdir    = "";

# --------------- #
# usage and input
# --------------- #
if ( $#ARGV < 2 ){
die("
usage: virusgathererHPC-3-assembly.pl <family> <project_id> <assembly_id> [options]\n
options:
\t-t=<integer>\tnumber of threads (default: $threads)
\t-b\t\tuse blastp for inital seed search (default: hmmsearch)
\t-v=<project>\tuse contigs produced by virushunter search as seed
\t\t\t(default: initial seed search by hmmsearch)
\t-se=<string>\tuser-provided seed file
\t-m=<integer>\trequested RAM in Mb (default: $memPerCPU)
\t-a=<string>\tGenSeed-HMM assembler parameter (default: $GS_assembler)
\t-u\t\tuse auxiliary starter assembler (cap3) for newbler (default: none)
\t-s=<integer>\tGenSeed-HMM seed size parameter (default: $GS_seed_size)
\t-e=<float>\tGenSeed-HMM blastn E-value parameter (default: $GS_blast_E)
\t-i=<integer>\tGenSeed-HMM blastn identity parameter (default: $GS_blast_ide)
\t-w=<integer>\tGenSeed-HMM blastn word size parameter (default: $GS_blast_ws)
\t-l=<integer>\tGenSeed-HMM max contig length parameter (default: $GS_contig_len)
\t-k=<integer>\tGenSeed-HMM max iterations parameters (default: $GS_max_itr)
\t-h=<integer>\tGenSeed-HMM alignment threshold parameter (default: $GS_aln_th)
\t-p=<integer>\tGenSeed-HMM CAP3 percent identity parameter (default: $cap_identity)
\t-ol=<integer>\tGenSeed-HMM CAP3 absolute overlap parameter (default: $cap_overlap)
\t-oh=<integer>\tGenSeed-HMM CAP3 percent overhang paraemter (default: $cap_overhang)
\t-di=<l|r|b>\tGenSeed-HMM expansion direction (default: $cap_overhang)
\t-TNT=<s|t>\tassemble full data set using SPAdes (s, under construction) or
\t\t\tTrinity (t, only recommended for RNA-seq data)
\t-lustre\t\tsave files on /lustre/ssd
\t-o=<string>\tuse output directory <o>
".
#\t-q=<yes|no>\tGenSeed-HMM use quality parameters (default: $GS_use_qual)
"\n");
}


# -------------- #
# get parameters
# -------------- #
my $family   = shift;
my $project  = shift;
my $assembly = shift;
while ( $#ARGV > -1 ){
	my $option = shift;
	if ( $option =~ /-t=(\d+)/ ){ $threads       = $1; }
	if ( $option eq "-b"       ){ $useBlastSeed  =  1; }
	if ( $option =~ /-v=(.+)/  ){ $useHitReads   = $1; }
	if ( $option =~ /-se=(.+)/ ){ $userSeed      = $1; }
	if ( $option =~ /-m=(\d+)/ ){ $memPerCPU     = $1; }
	if ( $option =~ /-a=(.+)/  ){ $GS_assembler  = $1; }
	if ( $option eq "-u"       ){ $GS_starter    =  1; }
	if ( $option =~ /-s=(\d+)/ ){ $GS_seed_size  = $1; }
	if ( $option =~ /-e=(.+)/  ){ $GS_blast_E    = $1; }
	if ( $option =~ /-i=(\d+)/ ){ $GS_blast_ide  = $1; }
	if ( $option =~ /-w=(\d+)/ ){ $GS_blast_ws   = $1; }
	if ( $option =~ /-l=(\d+)/ ){ $GS_contig_len = $1; }
	if ( $option =~ /-k=(\d+)/ ){ $GS_max_itr    = $1; }
	if ( $option =~ /-h=(\d+)/ ){ $GS_aln_th     = $1; }
	if ( $option =~ /-p=(\d+)/ ){ $cap_identity  = $1; }
	if ( $option =~ /-ol=(\d+)/){ $cap_overlap   = $1; }
	if ( $option =~ /-oh=(\d+)/){ $cap_overhang  = $1; }
	if ( $option =~ /-di=(.+)/ ){ $GS_direction  = $1; }
	if ( $option =~ /-TNT=(.+)/){ $TNT           = $1; $TNT = 1 if $TNT eq "t"; $TNT = 2 if $TNT eq "s"; }
	if ( $option eq "-lustre"  ){ $useLustre     =  1; }
	if ( $option =~ /-o=(.+)/  ){ $userOutdir    = $1; }
#	if ( $option =~ /-q=(.+)/  ){ $GS_use_qual   = $1; }
}


# -------------- #
# initialization
# -------------- #
my $resdir  = "$resdir0/$family/results/$project/assemblies";
   $tmpdir  = sprintf "%s/virusgatherer", virusHGutils::get_workspace_path( "lustre" )  if $useLustre == 1;
   $tmpdir .= "/$family/$project/$assembly";
   $tmpdir  = "$userOutdir/$assembly"  if $userOutdir ne "";

# log start time
my ($stime, $etime, $rtime, $ho, $mi);
$stime = time();


# ---------------------------------- #
# move data to working tmp directory
# ---------------------------------- #
my $orig_wdir = $tmpdir;
my $sraSize   = `du $tmpdir/data/*.sra | awk '{s+=\$1}END{print s}'`;  chomp( $sraSize );
   $sraSize   = $sraSize / 1024 / 1024;

my $tmpdirtmp = "";
   $tmpdirtmp = $tmpdir2  if ( $sraSize <= 25 );
   $tmpdirtmp = $tmpdir1  if ( $sraSize <= 10 );
if ( $tmpdirtmp ne "" ){
        my $diskfree = `df $tmpdirtmp | tail -n 1 | awk '{print \$4};'`;  chomp( $diskfree );
           $diskfree = $diskfree / 1024 / 1024;
        $tmpdirtmp = ""  if ( ($sraSize * $sizeInflatF) > $diskfree );
}

if ( $tmpdirtmp ne "" ){
        printf "[virusgatherer] moving data to temporary directory %s\n", $tmpdirtmp;
        `mkdir                   $tmpdirtmp/$family-$project-$assembly`;
        `cp -r       $tmpdir/*   $tmpdirtmp/$family-$project-$assembly/`;
        $tmpdir = sprintf "%s", "$tmpdirtmp/$family-$project-$assembly";
}


# -------------- #
# calculations
# -------------- #

# merge all fasta data files of the assembly
if ( ! -e "$tmpdir/data/data.fasta"  and  $TNT == 0 ){
	printf "[virusgatherer] merging fasta files for the assembly\n";
	`cat $tmpdir/data/*trim.q*.fasta > $tmpdir/data/data`;
}

# delete old data file
if ( ! -e "$tmpdir/data/data.fasta"  and  $TNT == 0 ){
	printf "[virusgatherer] deleting original data files\n";
	`mv $tmpdir/data/*.log $tmpdir/`;
	`rm $tmpdir/data/*.*`;
	`mv $tmpdir/data/data $tmpdir/data/data.fasta`;
}
if ( $TNT != 0 ){
	printf "[virusgatherer] deleting original data files\n";
	`rm $tmpdir/data/*_[1-2].fastq`;
	`rm $tmpdir/data/*.sra`;
}

# assemble
my $finalContigFile = "";
runGenseed()  if  $TNT == 0;
runTrinity()  if  $TNT == 1;

# blast final contigs against viral references
my $virusdbFile = "$virusdbdir/$family-db.fas";
if ( -e $finalContigFile ){
	my $blastcmd = "tblastx -query $finalContigFile -db $virusdbFile -num_threads $threads -max_hsps 1 -max_target_seqs 1 -outfmt '6 qseqid qlen evalue pident sseqid slen'";
	`$blastcmd > $finalContigFile-vs-virusDB-tblastx.tsv`;
}else{
	printf "Could not find final contig file '$finalContigFile'."
}


# -------- #
# clean up
# -------- #
if ( $tmpdirtmp ne "" ){
        printf "[virusgatherer] moving results to original working directory %s\n", $orig_wdir;
        #`cp -r  $tmpdir/results-* $orig_wdir/`;
        `rsync -avzr  $tmpdir/* $orig_wdir/`;
        `rm -rf $tmpdir`;
}



# log end time and calculate runtime
$etime = time();
$rtime = $etime - $stime;
$ho = int( $rtime / 3600 );  $rtime -= $ho * 3600;
$mi = int( $rtime /   60 );  $rtime -= $mi *   60;
printf "[virusgatherer] done - runtime: %d h %d min %d sec\n", $ho, $mi, $rtime;

#
# functions
# 

# run iterative assembly using GenSeed-HMM
sub runGenseed{
	# define seed for Genseed-HMM
	my $seedFile = "$profiledir/$family-seed.hmm";
	   $seedFile = "$querydir/$family-seed.fas" if ( $useBlastSeed == 1 );
	if ( $useHitReads ne "" ){
		printf "[virusgatherer] using virushunter contigs as seed\n";
		$seedFile = "$tmpdir/gen.seed.fas";
		my @datLogFiles = glob( "$tmpdir/*-data.txt" );
		my $hitDirsExist = 0;
		my $cwd = `pwd`;  chomp( $cwd );
		# concatenate seed contigs/singlets
		foreach my $dlf ( @datLogFiles ){
			$dlf =~ /.*\/([^\/]+)/;
			my @dlf = split( "-", $1 );
			my $hitdir  = "$resdir0/$family/results/$useHitReads/hits/".$dlf[0];
			my $hittar  = "$resdir0/$family/results/$useHitReads/hits.tar.gz";
			if ( ! -d $hitdir  and  -e $hittar ){
				chdir( $tmpdir );
				`tar -zxvf $hittar hits/$dlf[0]`;
				chdir( $cwd );
				$hitdir = "$tmpdir/hits/$dlf[0]/contigs.singlets.fas.gz";
				`zcat $hitdir >> $seedFile`;
				`rm -rf $tmpdir/hits`;
			}else{
				$hitdir .= "/contigs.singlets.fas.gz";
				next if ! -e $hitdir;
				`zcat $hitdir >> $seedFile`;
			}
			$hitDirsExist++;
		}
		# abort if no seed contigs/singlets found 
		if ( $hitDirsExist == 0 ){
			print "\nNo virushunter contigs found. Terminating.\n\n";
			exit();
		}
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
	}
	# run GenSeed-HMM
	printf "[virusgatherer] running GenSeed-HMM for iterative assembly\n";
	my $cmd  = "$GenSeedExec"; 
	   $cmd .= sprintf " -seed %s",				$seedFile;
	   $cmd .= sprintf " -db %s",				"$tmpdir/data/data.fasta";
	   $cmd .= sprintf " -assembler %s",			$GS_assembler;
	   $cmd .= sprintf " -ext_seed_size %d",		$GS_seed_size;
	   $cmd .= sprintf " -output %s",			"$tmpdir/results-$GS_assembler";
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
	   $cmd .= sprintf " -assembler_parameters \"-o $cap_overlap -h $cap_overhang -p $cap_identity\""  if ( $GS_assembler eq "cap3" );
	   $cmd .= sprintf " -aux_starter yes",     if ( $GS_starter == 1 );
	   $cmd .= sprintf " -exp_direction left",  if ( $GS_direction eq "l" );
	   $cmd .= sprintf " -exp_direction right", if ( $GS_direction eq "r" );
	#printf "[virusgatherer] command: %s\n", $cmd;
	`$cmd`;
	# final contig(s)
	$finalContigFile = "$tmpdir/results-$GS_assembler/final_result_dir/final_positive_contigs.fasta";
} # end runGenseed


# assemble full data set using Trinity
sub runTrinity{
	# prepare files and parameters
	my $maxMem  = int( $memPerCPU * $threads / 1024 );
	my @fqFiles = glob( "$tmpdir/data/*trim*.fastq" );
	my $cmd  = sprintf "trinity --seqType fq --CPU %d --max_memory %dG --output %s", $threads, $maxMem, "$tmpdir/results-trinity";
	   $cmd .= sprintf " --left %s --right %s", $fqFiles[0], $fqFiles[1]  if $#fqFiles == 1;
	   $cmd .= sprintf " --single %s",          $fqFiles[0],              if $#fqFiles == 0;
	# run Trinity
	printf "[virusgatherer] running Trinity\n";
	printf "$cmd\n";
	`$cmd > $tmpdir/trinity.log 2>&1`;
	# final contigs
	$finalContigFile = "$tmpdir/results-trinity/Trinity.fasta";
} # end runTrinity

