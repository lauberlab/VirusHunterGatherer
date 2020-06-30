#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

use lib "/home/lauber/lib/perl";
use Parallel::ForkManager;
use List::Util qw( shuffle );
use POSIX;

use lib "/home/lauber/lauber-2015-virusHunter";
use      virusHGutils;


# parameters
my $memPerCPU   = 1875; #2580;
my $cpusPerTask = 6; # deprecated
my $hitsPerRead = 2;
my $debugMode   = 0;
my $evalueB     = 1;
my $evalueH     = 10;
my $sizeLimit   = 50;
my $sizeTHlarge = 10;
my $sizeTHsmall = 1;
my $hoursPerGb  = 1.5; #0.5;
my $hoursRefSeq = 5;
my $hoursMaxHPC = 168; # submission of jobs running longer than 7 days are not allowed
my @ncbiIPs     = ("130.14.250.7","130.14.250.10","130.14.250.11",
                   "130.14.250.12","130.14.250.13");
my @ebiIPs      = ("fasp.sra.ebi.ac.uk");
my $EBIonly     = 0;
my $delayDL     = 0;
my $queryDir    = "/home/lauber/lauber-2015-virusHunter/queries";
my $profileDir  = "/home/lauber/lauber-2015-virusHunter/profiles";
my $motifDir    = "/home/lauber/lauber-2015-virusHunter/motifs";
my $doProfile   = 0;
my $doMotifs    = 0;
my $doFilter    = 0;
my $doRFblast   = 1;
my $filterHg38  = 0;
my $mapHg38     = 0;
my $maxSRAruns  = 10000000;
my $scLustre    = 6; # stripecount
my $scScratch   = 6; # stripecount
my $fastq_input = 0;
my $reservation = "";


# input
if ( $#ARGV < 3 ){
die("\n # # # # # # # # # # # # # # # # # # # # # # # # #
 # virushunter HPC version 2.30                  #
 #  identify viral sequences in SRA through      #
 #  sequence- or profile-based similarity search #
 # designer:  Stefan Seitz and Chris Lauber      #
 # developer: Chris Lauber                       #
 # email:     chris.lauber\@tu-dresden.de         #
 # # # # # # # # # # # # # # # # # # # # # # # # #\n
usage: virushunterHPC-submit.pl <family> <sra_accs_file> <project> <tasks> [options]\n
options:
"
#\t-t=<integer>\tuse <t> threads per task (default: all available)
."\t-e=<float>\tset E-value threshold to <e>
\t\t\t(default: $evalueB for Blast, $evalueH for HMMer)
\t-l=<float>\tlimit analysis to SRA run sizes of at most <l> Gb
\t\t\t(default: $sizeLimit)
\t-hr=<float>\tallocate a maximum of <hr> hours of runtime per Gb of data;
\t\t\tset to zero for no limit on runtime (default: $hoursPerGb)
\t-ha=<float>\tallocate an additional <ha> hours of runtime (for filtering)
\t\t\t(default: $hoursRefSeq)
\t-p\t\tdo profile-based search (default: Blast-based search)
\t-m\t\tperform additional motif search
\t-f\t\trun additional filtering step against known viral genomes
\t-d\t\tdebug mode (default: disabled)
\t-k=<integer>\tonly analyze first <k> SRA runs (default: all)
\t-dly=<float>\tdelay consecutive downloads for <dly> minutes if tasks==1
\t-rsb\t\tdo not run counter-blast search against RefSeq non-redundant proteins
\t-hg38f\t\tfilter against hg38 in filtering step 1
\t-hg38m\t\tmap reads against hg38 genes (meaningful only for transcriptomes)
\t\t\t(only recommended for analyses of human SRA experiments)
\t-fastq\t\tprovided are not SRA identifiers for download but names of
\t\t\tlocal fastq.gz files
\t-ebi\t\tdownload only from EBI (default: from NCBI and EBI)
\t-rsv=<string>\tuse reservation <rsv> on Taurus (default: submit to normal queue)
\n");
}

my $family      = shift;
my $sraFile     = shift;
my $project     = shift;
my $cores       = shift;
while ( $#ARGV >= 0 ){
	my $option = shift;
	if ( $option =~ /-t=(\d+)/     ){  $cpusPerTask = $1; }
	if ( $option =~ /-e=(.*)/      ){  $evalueB     = $1; $evalueH = $1; }
	if ( $option =~ /-l=(\d+.*)/   ){  $sizeLimit   = $1; }
	if ( $option =~ /-hr=(\d+.*)/  ){  $hoursPerGb  = $1; }
	if ( $option =~ /-ha=(\d+.*)/  ){  $hoursRefSeq = $1; }
	if ( $option =~ /-dly=(\d+.*)/ ){  $delayDL     = $1; }
	if ( $option =~ /-k=(\d+)/     ){  $maxSRAruns  = $1; }
	if ( $option eq "-p"           ){  $doProfile   =  1; }
	if ( $option eq "-m"           ){  $doMotifs    =  1; }
	if ( $option eq "-f"           ){  $doFilter    =  1; }
	if ( $option eq "-d"           ){  $debugMode   =  1; }
	if ( $option eq "-rsb"         ){  $doRFblast   =  0; }
	if ( $option eq "-hg38f"       ){  $filterHg38  =  1; }
	if ( $option eq "-hg38m"       ){  $mapHg38     =  1; }
	if ( $option eq "-fastq"       ){  $fastq_input =  1; }
	if ( $option eq "-ebi"         ){  $EBIonly     =  1; }
	if ( $option =~ /-rsv=(.+)/    ){  $reservation = $1; }
}


# global parameters
#my $bindir     = "/home/lauber/lauber-2015-virusHunter/developer";
my $bindir     = "/home/lauber/lauber-2015-virusHunter";
my $tmpdir     = sprintf "%s/lauber-2015-virusHunter/$family",  virusHGutils::get_workspace_path( "scratch" );
my $logtmpdir  = sprintf "%s/virushunter/$family/$project",     virusHGutils::get_workspace_path( "lustre"  );
my $resdir     = "/projects/p_sra/$family/results";
#my $logdir     = "/projects/p_sra/logs";
my $logdir     = "/projects/p_sra/$family/results/$project/log";
my $statdir    = "/projects/p_sra/$family/statistics";
my %jobs       = ();

# current directory
my $cmddir = `pwd`;  chomp( $cmddir );

# create family subdirectories if not alread present
if ( ! -d "$tmpdir" ){
	mkdir( "$tmpdir" );
	`lfs setstripe -c 6 $tmpdir`; # set stripecount to balance load between OSTs (due to large file sizes) 
}
if ( ! -d "$tmpdir/$project" ){
	mkdir( "$tmpdir/$project" );
}
if ( ! -d "/projects/p_sra/$family" ){
	mkdir(  "/projects/p_sra/$family" );
	mkdir(  "/projects/p_sra/$family/results" );
	mkdir(  "/projects/p_sra/$family/statistics" );
	open(C,">/projects/p_sra/$family/statistics/data.txt"); close(C);
	open(C,">/projects/p_sra/$family/statistics/search.txt"); close(C);
}
if ( ! -d sprintf "%s/virushunter/$family", virusHGutils::get_workspace_path( "lustre" ) ){
	my $stripe_dir = sprintf "%s/virushunter/$family", virusHGutils::get_workspace_path( "lustre" );
	mkdir( $stripe_dir );
	`lfs setstripe -c 6 $stripe_dir`; # set stripecount to balance load between OSTs (due to large file sizes) 
}
if ( ! -d "$logtmpdir" ){
	mkdir(  "$logtmpdir" );
}
if ( ! -d "$logtmpdir/log" ){
	mkdir(  "$logtmpdir/log" );
}
$logtmpdir .= "/log";

# create project dir if not already present
if ( ! -d "$resdir/$project" ){
	# folders
	mkdir( "$resdir/$project" );
	mkdir( "$resdir/$project/hits" );
	mkdir( "$resdir/$project/log" );
	# empty file with SRA identifiers that were already processed
	open(C,">$resdir/$project/log/completed.txt"); close(C);
	# empty archive for log files
	`tar -cf $resdir/$project/log/log.tar -T /dev/null;  gzip $resdir/$project/log/log.tar`;
	# create zip archives
	chdir( "$resdir/$project" );
	`tar -rf hits.tar hits`;  `gzip hits.tar`;  `rm -rf hits`;
	`tar -rf  log.tar  log`;  `gzip  log.tar`;  `rm -rf  log`;
}

# unpack completed
chdir( "$resdir/$project" );
`rm  -rf log`;
`tar -zxvf log.tar.gz log/completed.txt`;
chdir( $cmddir );


# determine query
my $querySeqs   = "$queryDir/$family-query.fas";
my $profile     = "$profileDir/$family-profile.hmm";
my $motifs      = "$motifDir/$family-motifs.meme";
my $filterSeqs  = "$queryDir/$family-filter.fas";

# verify that query/profile file exists
if ( $doProfile == 0 and ! -e $querySeqs ){
	die( "\nNo query sequence file for family '$family' found!\n\n" );
}
if ( $doProfile == 1 and ! -e $profile ){
	die( "\nNo query profile file for family '$family' found!\n\n" );
}

# read SRA IDS from file and determine its size
# do not run SRAs that have already been processed previously
# do not run SRAs that exceed the size limit
my %srasall   = ();
my $fastq_dir = "";
open( SRA,"<$sraFile");
while ( my $line = <SRA> ){
	# SRA run id
	chomp($line);  next if $line eq "";
	if ( $fastq_input ){
		if ( $line =~ /.*\/.*/ ){
			$line =~ /(.*)\/([^\/]+)\.fastq\.gz/;
			$fastq_dir = $1;
			$line      = $2;
		}else{
			$fastq_dir = ".";
		}
	}
	# SRA run already processed ?
	#next if -e "$resdir/$project/log/$line-tblastn.out";
	#next if -e "$resdir/$project/log/$line-hmmsearch.out";
	#my $findB = `if tar -tf $resdir/$project/log/log.tar.gz $line-tblastn.out   > /dev/null 2>&1; then echo 1; else echo 0; fi`;  chomp( $findB );
	#my $findH = `if tar -tf $resdir/$project/log/log.tar.gz $line-hmmsearch.out > /dev/null 2>&1; then echo 1; else echo 0; fi`;  chomp( $findH );
	#next if ( $findB or $findH );
	my $findC = `if grep '$line' $resdir/$project/log/completed.txt > /dev/null 2>&1; then echo 1; else echo 0; fi`;  chomp( $findC );
	next if $findC;
	# include this run
	$srasall{ $line } = 1;
	last if ( scalar keys %srasall >= $maxSRAruns );
}
close(SRA);

# get size of each SRA
my $sras = {};
my ( $sizeMin, $sizeMean, $sizeMax, $sizeRm ) = ( 1e6, 0, 0, 0 );
my ( $wget_cmd, $wget_res, $wget_n, $wget_fail ) = ( "", "", 0, 0 );
my @srasHalf1 = ();  my @srasHalf2 = ();
if ( $fastq_input == 0 )
{
my @srasall = keys %srasall;
printf "%s retrieving metadata for %d SRA runs\n", printTime(), scalar @srasall;
while ( $#srasall > -1 ){
	my $sraid = pop( @srasall );
	# SRA run above size limit ? (100 runs in a batch)
	if ( scalar keys %srasall == 1 ){
		$wget_cmd  = "wget -q -O - 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=$sraid";
	}
	#if ( $wget_n == 100 or $#srasall < 0 ){
	if ( $wget_n ==  50 or $#srasall < 0 ){
		$wget_cmd .= " OR $sraid" if ( $#srasall < 0 and scalar keys %srasall > 1 );
		$wget_cmd .= "'";
		# get info for this batch
		my $wget_itr = 0;
		   $wget_res = "";
		   $wget_n   = 0;
		while ( $wget_res eq "" ){
			#printf "\n%s\n\n", $wget_cmd;
			$wget_res = `$wget_cmd`;
			$wget_itr++;
			last if $wget_itr > 5;
			sleep( 12 ) if $wget_res eq "";
		}
		if ( $wget_res ne "" ){
			my @lin = split /\n/, $wget_res;
			my @hdr = split /,/,  $lin[0];
			for ( my $si=1; $si<=$#lin; $si++ ){
				my $srasize = 2048; # init to about-average size of 2 Gb
				my $sid     = "";
				my @val = split /,/,  $lin[ $si ];
				for ( my $wget_i=0; $wget_i<=$#hdr; $wget_i++ ){
					$srasize = $val[ $wget_i ] if ( $hdr[ $wget_i ] eq "size_MB" );
					$sid     = $val[ $wget_i ] if ( $hdr[ $wget_i ] eq "Run" );
				}
				next if $srasize !~ /^[0-9,.E]+$/;
				$srasize /= 1024;
				next if $sid eq "";
				next if exists $sras->{$sid};
				next if ! exists $srasall{ $sid };
				if ( $srasize > $sizeLimit){
					$sizeRm++;
					next;
				}
				# some size statistics
				$sizeMin   = $srasize if $srasize < $sizeMin;
				$sizeMax   = $srasize if $srasize > $sizeMax;
				$sizeMean += $srasize;
				# run this SRA run
				$sras->{$sid} = {};
				$sras->{$sid}->{'size'} = $srasize;
				# put in first or second analysis batch depending on SRA size
				if ( $srasize >= $sizeTHlarge or $srasize < $sizeTHsmall ){  push( @srasHalf1, $sid );
				}else{                                                       push( @srasHalf2, $sid ); }
			}
		}else{
			$wget_fail += $wget_n;
		}
		# init next batch
		$wget_cmd  = "wget -q -O - 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=$sraid";
	}elsif ($wget_cmd eq ""){
		$wget_cmd  = "wget -q -O - 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=$sraid";
	}else{
		$wget_cmd .= " OR $sraid";
	}
	$wget_n++;
}
printf "%s %d SRA runs exceed the size limit of %.1f Gb\n", printTime(), $sizeRm, $sizeLimit;
printf "%s %d SRA failed annotation retrieval from SRA\n",  printTime(), $wget_fail  if ( $wget_fail > 0 );
}

# order SRAs taking into account size
#  randomly mixing very small and very large SRAs in first half
#  and remaining SRAs in second half
my @sraOrdered = ();
if ( $fastq_input ){
	@sraOrdered = keys %srasall;
}else{
	if ( $sizeMax < $sizeTHlarge ){
		@sraOrdered = ( @srasHalf1, @srasHalf2 );
		@sraOrdered = shuffle @sraOrdered;
	}else{
		@srasHalf1  = shuffle @srasHalf1;
		@srasHalf2  = shuffle @srasHalf2;
		@sraOrdered = ( @srasHalf1, @srasHalf2 );
	}
}

# set provider and IP for data download
my ($ipiNCBI, $ipiEBI)  = ( -1, -1 );
my $provider = "ncbi";
   $provider = "ebi"  if ( $EBIonly );
my $sraorder = 0;
foreach my $sraid ( @sraOrdered ){
	## provider for data download
	#if ( $provider eq "" ){
	#	if ( $useNCBI == 1 ){  $provider = "ncbi"; }
	#	else{                  $provider = "ebi";  }
	#}else{
	#	if    ( $provider eq "ncbi"  and  $useEBI  == 1 ){  $provider = "ebi";  }
	#	elsif ( $provider eq "ebi"   and  $useNCBI == 1 ){  $provider = "ncbi"; }
	#}
	$sras->{$sraid}->{'provider'} = $provider;
	# IP/hostname for data download
	if ( $provider eq "ncbi" ){
		$ipiNCBI++;
		$ipiNCBI = 0 if $ipiNCBI > $#ncbiIPs;
		$sras->{$sraid}->{'ip'} = $ncbiIPs[ $ipiNCBI ];
	}else{
		$ipiEBI++;
		$ipiEBI = 0 if $ipiEBI > $#ebiIPs;
		$sras->{$sraid}->{'ip'} = $ebiIPs[ $ipiEBI ];
	}
	# remember order
	$sraorder++;
	$sras->{$sraid}->{'ni'} = $sraorder;
}

# function to report current time
sub printTime{
	my @date = split / +/, `date`;
	return sprintf "[%s-%s %s] ", $date[1], $date[2], $date[3];
}

# print progress header
if ( $fastq_input == 0 ){
	my $prgs = 0;
	if ( scalar keys %{$sras} == 0 ){  print  "no new SRA runs to be processed\n";  exit; }
	my $sizeTotal = sprintf "%.1f Gb", $sizeMean;
	   $sizeTotal = sprintf "%.1f Tb", $sizeMean / 1024 if $sizeMean >= 1024;
	printf "%s processing %d SRA runs comprising %s in total (min: %1.f Gb | mean: %.1f Gb | max: %.1f Gb)\n", 
		printTime(), scalar keys %{$sras}, $sizeTotal, $sizeMin, $sizeMean / scalar keys %{$sras}, $sizeMax;
	printf "                             |";
	while ( 1 ){
		$prgs += 250;
		last if $prgs > scalar keys %{$sras};
		if ( $prgs % 1000 == 0 ){  print "|";
		}else{                     print "_"; }
	}
	print "\n                   progress:  ";
}

# copy completed, data, and search logfiles to tmpdir
# (because of flock performance) and make backups
`cp $resdir/$project/log/completed.txt $resdir/$project/log/completed.txt.backup`;
`cp $resdir/$project/log/completed.txt $tmpdir/virushunterHPC-$project-completed.txt`;
`cp $statdir/data.txt   $statdir/data.txt.backup`;
`cp $statdir/data.txt   $tmpdir/virushunterHPC-$project-data.txt`;
`cp $statdir/search.txt $statdir/search.txt.backup`;
`cp $statdir/search.txt $tmpdir/virushunterHPC-$project-search.txt`;

# temporary file to store job IDs
my $jobidfile = "$tmpdir/virushunterHPC-$project-jobIDs.txt";
open(J,">$jobidfile"); close(J);

# fork into <cores> child processes, download data in parallel
# and submit search job for each process
# init fork
my $pm = new Parallel::ForkManager($cores);
foreach my $sid ( @sraOrdered ){
	# start fork
	my $pid = $pm->start and next;
	# download or copy data
	my $cmd = "";
	if ( $fastq_input ){
		   $cmd  = sprintf "cp %s/%s %s/", $fastq_dir, $sid.".fastq.gz", $tmpdir;
	}else{
		my $ip   = $sras->{$sid}->{'ip'};
		my $prov = $sras->{$sid}->{'provider'};
		   $cmd  = sprintf "%s/virushunterHPC-1-data.pl %s %s %s %s %s -l=%s", $bindir, $family, $project, $sid, $ip, $prov, $sizeLimit;
		   $cmd .= " -p"             if ( $doProfile == 1 );
		   $cmd .= " -dly=$delayDL"  if ( $cores == 1 );
	}
	`$cmd > $logtmpdir/$sid-data.out 2> $logtmpdir/$sid-data.err`;
	# determine size of run and set max runtime accordingly
	my $maxRT = 0;
	if ( -e "$tmpdir/$sid.sra" or -e "$tmpdir/$sid.fastq.gz" ){
		my $sraSize  = 0;
		if ( $fastq_input ){
			$sraSize  = `ls -l $tmpdir/$sid.fastq.gz | awk '{print \$5}'`;
		}else{
			$sraSize  = `ls -l $tmpdir/$sid.sra      | awk '{print \$5}'`;
		}
		   $sraSize /= ( 1024 * 1024 * 1024 );
		   $maxRT    = $hoursPerGb;  $maxRT *= 2 if $doProfile == 1;
		   $maxRT    = ceil( $maxRT * $sraSize + $hoursRefSeq );
		   $maxRT    = $hoursMaxHPC  if ( $maxRT > $hoursMaxHPC );
	}
	# submit search job
	my $jobid = submitJob( $sid, $maxRT );
	if ( $jobid eq "" ){
		printf STDERR "unable to submit job for '$sid'\n";
	}else{
		open( JID, ">>$jobidfile") or die("$0: could not open file: $!\n");
		flock(JID,2) or die("$0: could not lock file: $!\n");
		print JID "$jobid\n";
		close(JID) or die("$0: could not close file: $!\n");
	}
	# report some progress
	print "." if ( $sras->{$sid}->{'ni'} % 250 == 0 );
	# end fork
	$pm->finish;
}
# wait for all forks
$pm->wait_all_children;
print "\n";

# do not wait for completion of all jobs anymore
# it means that 'manual cleanup' has to be run afterwards
printf "%s will not wait for jobs to finish ('manual cleanup' required)\n", printTime();
printf "%s done\n", printTime();
exit;


# wait for all jobs to finish
printf "%s waiting for all jobs to finish\n", printTime();
my %jobids = ();
open( JIDS, "<$jobidfile" ) or die("$0: could not open file: $!\n");
while (my $line = <JIDS> ){
	chomp($line);  next if $line eq "";
	$jobids{$line} = 1;
}
close(JIDS) or die("$0: could not close file: $!\n");
#while (1){
#	last if ( scalar keys %jobids == 0 );
#	my $squeue = `squeue -u lauber -h -o %i -t R,PD 2>&1`;
#	if ( $squeue =~ /.*error.*/ ){
#		`echo "$squeue" > $tmpdir/virushunterHPC-$project-squeue.err`;
#	}else{
#		my @squeue = split( /\n/, $squeue );
#		my %squeue = ();
#		foreach ( @squeue ){ $squeue{$_} = 1; }
#		foreach ( keys %jobids ){
#			#my $squeue = `squeue -u lauber 2> /dev/null | grep $_ 2> /dev/null`;
#			#delete $jobids{$_} if $squeue eq "";
#			next if exists $squeue{ $_ };
#			delete $jobids{ $_ };
#		}
#	}
#	sleep(30);
#}

open( SQUEUE, ">$tmpdir/virushunterHPC-$project-squeue.out" )  if ($debugMode != 0);
while (1){
	sleep(60);
	my ($running, $unknown) = (0,0);
	foreach ( keys %jobids ){
		my $sacct = "";
		   $sacct = `sacct -b -n -j $_ 2>&1`;
		printf SQUEUE "------\nsacct:\n------\n$sacct"  if ($debugMode != 0);
		if ( $sacct =~ /.*error.*/  or  $sacct eq "" ){
			`echo "$sacct" > $tmpdir/virushunterHPC-$project-squeue.err`;
			$unknown++;
			sleep( 60 );
		}else{
			my @sacct = split( /\s+/, $sacct );
			$jobids{$_} = $sacct[1];
			if ( $sacct[1] eq "RUNNING" or $sacct[1] eq "PENDING" ){
				$running++;
			}else{
				printf SQUEUE "%s finished\n", $_  if ($debugMode != 0);
			}
		}
	}
	last if ($running == 0  and  $unknown == 0);
}
close(SQUEUE)  if ($debugMode != 0);

# copy results from temporary to final project directory
`cp -r $tmpdir/$project/* $resdir/$project/hits/`;
if ( $debugMode == 0 ){  `rm -rf $tmpdir/$project`; }

# copy log-files from temporary to final project directory
#`cp $logtmpdir/* $logdir/`;
chdir( $logtmpdir );
`ls * > logfiles2copy.txt`;
`gunzip  $logdir/log.tar.gz`;
`tar -rf $logdir/log.tar -T logfiles2copy.txt`;
`gzip    $logdir/log.tar`;
chdir( $bindir );
if ( $debugMode == 0 ){  `rm -rf $logtmpdir`; }

# copy back completed, data, and search stat-files
`cp $tmpdir/virushunterHPC-$project-completed.txt $resdir/$project/log/completed.txt`;
`cp $tmpdir/virushunterHPC-$project-data.txt      $statdir/data.txt`;
`cp $tmpdir/virushunterHPC-$project-search.txt    $statdir/search.txt`;

# DONE
printf "%s done\n", printTime();

# function to compile job file for the search part and submit it
sub submitJob {
	my $sid = shift;
	my $mrt = shift;
	# which machines to use
	my $machines = "haswell,sandy,romeo";
	   $machines = "haswell" if ( $reservation ne "" );	
	#  compile job file for the search
	my $slurm = "#!/bin/bash\n";
	$slurm .= sprintf "#SBATCH --output=%s\n",		$logtmpdir."/$sid-search.out";
	$slurm .= sprintf "#SBATCH --error=%s\n",		$logtmpdir."/$sid-search.err";
	if ( $mrt > 0 )
	{
	$slurm .= sprintf "#SBATCH --time=%d:00:00\n",		$mrt; # in hours
	}
	$slurm .= sprintf "#SBATCH --nodes=%d\n",		1;
	$slurm .= sprintf "#SBATCH --ntasks-per-node=%d\n",	1;
	#$slurm .= sprintf "#SBATCH --cpus-per-task=%d\n",	$cpusPerTask;
	#$slurm .= sprintf "#SBATCH -N 1 --exclusive\n";
	$slurm .= sprintf "#SBATCH --exclusive=lauber\n";
	$slurm .= sprintf "#SBATCH --mem-per-cpu=%d\n",		$memPerCPU;
	$slurm .= sprintf "#SBATCH -p %s\n",			$machines;
	#$slurm .= sprintf "#SBATCH -L scratch\n";
        if ( $reservation ne "" ){
        $slurm .= sprintf "#SBATCH --reservation=%s\n",         $reservation;
        }
	$slurm .= sprintf "#SBATCH --job-name=%s\n\n",		"h".$sid."-search";
	$slurm .= sprintf "module load EMBOSS\n";
	#$slurm .= sprintf "module load modenv/both\n";
	#$slurm .= sprintf "module load mysql/6.0.11-client\n";
	$slurm .= sprintf "cd %s\n", $bindir;
	$slurm .= "srun --cpu_bind=none ";
	# parameters
	if ( $doProfile == 1 )
	{
	$slurm .= sprintf "./virushunterHPC-2-search.pl %s %s %s %s -t=\$SLURM_JOB_CPUS_PER_NODE -h=%d -e=%s -p",
                          $family, $profile,   $project, $sid, $hitsPerRead, $evalueH;
	}else
	{
	$slurm .= sprintf "./virushunterHPC-2-search.pl %s %s %s %s -t=\$SLURM_JOB_CPUS_PER_NODE -h=%d -e=%s",
                          $family, $querySeqs, $project, $sid, $hitsPerRead, $evalueB;
	}
	# optional parameters
	$slurm .= " -d"             if $debugMode   == 1;
	$slurm .= " -m=$motifs"     if $doMotifs    == 1;
	$slurm .= " -f=$filterSeqs" if $doFilter    == 1;
	$slurm .= " -rsb"           if $doRFblast   == 0;
	$slurm .= " -fastq"         if $fastq_input == 1;
	$slurm .= " -hg38f"         if $filterHg38  == 1;
	$slurm .= " -hg38m"         if $mapHg38     == 1;
	# create job file
	my $jobfile = "$tmpdir/$sid-search.slurm";
	open (SLURM,">$jobfile");
	print SLURM $slurm."\n";
	close(SLURM);
	# submit command
	my $cmd = "sbatch $jobfile 2> /dev/null | awk '{print \$4}'";
	my $jid = "";
	for ( my $try=0; $try<10; $try++ ){
		# submit job file and log job ID
		$jid = `$cmd 2> /dev/null`;  chomp($jid);
		if ( $jid eq "" ){  sleep( 60 );
		}else{		    last; }
	}
	# delete job file
	`rm -rf $jobfile` if $jid ne "";
	# return job ID
	return( $jid );
}

