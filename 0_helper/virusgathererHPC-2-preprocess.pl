#!/usr/bin/perl

##########################################################
# extract SRA files to fastq, followed by adapter
#  cutting and quality trimming
#
# author: Chris Lauber
##########################################################

# ------------ #
# load modules
# ------------ #
use warnings;
use strict;
use Cwd;
use Data::Dumper;

use lib "/home/lauber/lauber-2015-virusHunter";
use      virusHGutils;


# --------- #
# parameter
# --------- #
my $readLenTH   = 10;			# remove (technical) fastq files with too short reads
my @adaps3p     = ("AGATCGGAAGAG",	# Illumina Universal Adapter
                   "ATGGAATTCTCG",	# Illumina Small RNA Adapter
                   "CTGTCTCTTATA",	# Nextera Transposase Sequence
  	  	   "CGCCTTGGCCGT",	# SOLID Small RNA Adapter
		   "CGACCTGAGACTGCCAAGGCACACAGGGGATAGG"
		   );
my @adaps5p	= ("AAGCAGTGGTATCAACGCAGAGTGGCCATTACGGCCGGG", 			# 454
		   "ACACGTAGTATAAGCAGTGGTATCAACGCAGAGTTTTTGTTTTTTTCTTTTTTTTTT",
		   "ACACGTAGTATAAGCAGTGGTATCAACGCAGAGTTTTTGTTTTTCTTTTTT",
		   "ACACGTAGTATAAGCAGTGGTATCAACGCAGAGT",
		   "ACACGTAGTATAAGCAGTGGTATCAACGCAGAGGTACGCGGGG",
		   "AGTGGTATCAACGCAGAGTACGCGGGG",
		   "ACACGTAGTAT",
		   "TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG"			# 454 titanium
	           );
my $readQcut    = 20; #25;
my $readMinL    = 20;
my $threads     = 1;
my $TNT         = 0; # TNT mode: assembly of full data set
my $fastq_input = 0;
my $tmpdir1     = "/tmp";
my $tmpdir2     = sprintf "%s/virusgatherer", virusHGutils::get_workspace_path( "lustre" );
my $sizeInflatF = 7; # determined empirically


# --------------- #
# usage and input
# --------------- #
if ( $#ARGV < 0 ){
die("
usage: virusgathererHPC-2-preprocess.pl <sra_file> [options]\n
options:
\t-q=<integer>\ttrim 3' read ends with quality lower <q> (default: $readQcut)
\t-m=<integer>\tminium read length after trimming (default: $readMinL)
\t-t=<integer>\tnumber of threads (default: $threads)
\t-tnt\t\tTNT mode: no fastq-to-fasta transformation
\t-fastq\t\tinput is not SRA runs but local fastq.gz files
\n");
}


# -------------- #
# get parameters
# -------------- #
my $sraFile = shift;
while ( $#ARGV > -1 ){
	my $option = shift;
	if ( $option =~ /-q=(\d+)/ ){ $readQcut    = $1; }
	if ( $option =~ /-m=(\d+)/ ){ $readMinL    = $1; }
	if ( $option =~ /-t=(\d+)/ ){ $threads     = $1; }
	if ( $option eq "-tnt"     ){ $TNT         =  1; }
	if ( $option eq "-fastq"   ){ $fastq_input =  1; }
}


# ---------------------------------- #
# move data to working tmp directory
# ---------------------------------- #
$sraFile =~ /(.*)\/([^\/]+)$/;
my $orig_wdir = $1;
my $fileid    = $2;
my $sraid     = $fileid;  $sraid =~ s/\..*//;
my $sraSize   = `du $sraFile | awk '{print \$1}'`;  chomp( $sraSize );
   $sraSize   = $sraSize / 1024 / 1024;

my $tmpdir = "";
   $tmpdir = $tmpdir2  if ( $sraSize <= 25 );
   $tmpdir = $tmpdir1  if ( $sraSize <= 10 );
if ( $tmpdir ne "" ){
	my $diskfree = `df $tmpdir | tail -n 1 | awk '{print \$4};'`;  chomp( $diskfree );
	   $diskfree = $diskfree / 1024 / 1024;
	$tmpdir = ""  if ( ($sraSize * $sizeInflatF) > $diskfree );
}

if ( $tmpdir ne "" ){
	printf "[virusgatherer] moving data to temporary directory %s\n", $tmpdir;
	`mkdir          $tmpdir/$sraid`;
	`cp    $sraFile $tmpdir/$sraid/`;
	$sraFile = sprintf "%s/%s", "$tmpdir/$sraid/", $fileid;
}


# -------------- #
# calculations
# -------------- #

# log start time
my ($stime, $etime, $rtime, $ho, $mi);
$stime = time();

# get directory from file name
my ($pdir, $rid) = ("", "");
if ( $fastq_input ){
	$sraFile =~ /(.*)\/([^\/]+)$/;
	$pdir = $1;  $rid = $2;
	foreach my $fastqgz ( glob( "$sraFile/*.fastq.gz" ) ){
		if ( $TNT == 1 ){
			$fastqgz =~ s/\.gz$//;
			`gunzip $fastqgz.gz`;
			`mv $fastqgz $pdir/`;
		}else{
			my $zcat_cmd = sprintf "zcat %s >> %s/%s_1.fastq", $fastqgz, $pdir, $rid;
			`$zcat_cmd`;
		}
	}
	#`rm -rf $sraFile`;
}else{
	$sraFile =~ /(.*)\/([^\/]+)\.sra$/;
	$pdir = $1;  $rid = $2;
	printf "[virusgatherer] processing run $rid in directory $pdir\n";
	# extract SRA to fastq
	printf "[virusgatherer] extracting fastq from SRA run\n";
	my $cmd = "fastq-dump --skip-technical --split-files --readids -O $pdir -B $sraFile";
	   $cmd = "fastq-dump --skip-technical --split-files --defline-seq '\@\$ac.\$sn/\$ri' --defline-qual '+\$ac.\$sn/\$ri' -O $pdir -B $sraFile"  if $TNT == 1;
	#print "$cmd\n";
	`$cmd`;
}

# get read files
my @readFs = ();
foreach ( glob( "$pdir/$rid"."_*.fastq" ) ){
	if ( $_ =~ /^$pdir\/$rid\_\d+\.fastq$/ ){
		push( @readFs,$_ );
	}
}

# remove technical / too short reads
printf "[virusgatherer] checking read lengths for each read file\n";
foreach my $fq ( @readFs ){
	# get largest length of the first 1000 reads
	my ($ri, $maxl) = (0, 0);
	open( RF, "<$fq" );
	while ( my $line = <RF> ){
		chomp($line);  next if $line eq "";
		$ri++;
		if ( ($ri-2) % 4 == 0 ){
			if ( length($line) > $maxl ){  $maxl = length($line); }
		}
		last if $ri >= 4000;
	}
	close(RF);
	# remove if below size threshold
	if ( $maxl < $readLenTH ){
		`rm $fq`;
		$fq =~ /.*\/(.*)/;  $fq = $1;
		printf "[virusgatherer] removing %s, as reads are shorter than %d nt (likely technical reads)\n", $fq, $readLenTH;
	}
}

# update read files
@readFs = ();
foreach ( glob( "$pdir/$rid"."_*.fastq" ) ){
	if ( $_ =~ /^$pdir\/$rid\_\d+\.fastq$/ ){
		push( @readFs,$_ );
	}
}


# cut adapters and trim low quality bases (cutadapt)
#printf "[virusgatherer] cutting adapter sequences and trimming low-quality 3' ends for run '%s'\n", $rid;
#my $cut_log = $pdir."/$rid-cutadapt.log";
# paired or single reads?
#my @readNs = ();
#foreach ( @readFs ){
#	open(WC,"wc -l $_ |"); my $wcl = <WC>;  close(WC);
#	$wcl =~ /(\d+)\s.*/;
#	push( @readNs, $1 );
#}
#if ( $#readFs == 1 and ($readNs[0]-$readNs[1]) == 0 ){
#	printf "[virusgatherer] assuming paired reads\n", $rid;
#	my $cmd = "cutadapt";
#	#foreach ( @adaps ){  $cmd .= " -b $_"; }					
#	#foreach ( @adaps ){  $cmd .= " -B $_"; }
#	foreach ( @adaps3p ){  $cmd .= " -a $_"; }					
#	foreach ( @adaps3p ){  $cmd .= " -A $_"; }
#	foreach ( @adaps5p ){  $cmd .= " -g $_"; }					
#	foreach ( @adaps5p ){  $cmd .= " -G $_"; }
#	my $ofile = $readFs[0];  $ofile =~ s/\.fastq/\.trim\.q$read3q\.fastq/;
#	my $pfile = $readFs[1];  $pfile =~ s/\.fastq/\.trim\.q$read3q\.fastq/;
#	$cmd .= " -o $ofile";
#	$cmd .= " -p $pfile";		
#	#$cmd .= sprintf " -n %d -m %d -q %d", $#adaps+1, $readMinL, $read3q;
#	$cmd .= sprintf " -m %d -q %d", $readMinL, $read3q;
#	foreach ( @readFs ){  $cmd .= " $_"; }
#	`$cmd > $cut_log`;
#}else{
#	printf "[virusgatherer] assuming single reads\n", $rid;
#	for ( my $i=0; $i<=$#readFs; $i++ ){
#		my $cmd = "cutadapt";
#		#foreach ( @adaps ){  $cmd .= " -b $_"; }					
#		foreach ( @adaps3p ){  $cmd .= " -a $_"; }					
#		foreach ( @adaps5p ){  $cmd .= " -g $_"; }					
#		my $ofile = $readFs[$i];  $ofile =~ s/\.fastq/\.trim\.q$read3q\.fastq/;
#		$cmd .= " -o $ofile";
#		#$cmd .= sprintf " -n %d -m %d -q %d", $#adaps+1, $readMinL, $read3q;
#		$cmd .= sprintf " -m %d -q %d", $readMinL, $read3q;
#		$cmd .= " $readFs[$i]";
#		if ( $i == 0 ){  `$cmd  > $cut_log`;
#		}else{           `$cmd >> $cut_log`; }
#	}
#}				

# cut adapters and trim low quality bases (autoadapt)
printf "[virusgatherer] cutting adapter sequences and trimming low-quality 3' ends for run '%s'\n", $rid;
# paired or single reads ?
my @readNs = ();
foreach ( @readFs ){
	open(WC,"wc -l $_ |"); my $wcl = <WC>;  close(WC);
	$wcl =~ /(\d+)\s.*/;
	push( @readNs, $1 );
}
my $cut_log = $pdir."/$rid-autoadapt.log";
if ( $#readFs == 1 and ($readNs[0]-$readNs[1]) == 0 ){
	printf "[virusgatherer] assuming paired reads\n", $rid;
	my $infile1  = $readFs[0];
	my $infile2  = $readFs[1];
	my $outfile1 = $readFs[0];  $outfile1 =~ s/\.fastq/\.trim\.q$readQcut\.fastq/;
	my $outfile2 = $readFs[1];  $outfile2 =~ s/\.fastq/\.trim\.q$readQcut\.fastq/;
	my $cut_cmd = "autoadapt --threads=$threads --quality-cutoff=$readQcut --minimum-length=$readMinL";
	$cut_cmd .= " $infile1 $outfile1 $infile2 $outfile2";
	my $wdir = cwd;
	chdir( $pdir );
	`$cut_cmd > $cut_log`;
	chdir( $wdir );
}else{
	printf "[virusgatherer] assuming single reads\n", $rid;
	for ( my $i=0; $i<=$#readFs; $i++ ){
		my $infile  = $readFs[$i];
		my $outfile = $readFs[$i];  $outfile =~ s/\.fastq/\.trim\.q$readQcut\.fastq/;
		my $cut_cmd = "autoadapt --threads=$threads --quality-cutoff=$readQcut --minimum-length=$readMinL";
		$cut_cmd .= " $infile $outfile";
		my $wdir = cwd;
		chdir( $pdir );

		`$cut_cmd > $cut_log`;
		chdir( $wdir );
	}
}

# transform to fasta
if ( $TNT == 0 ){
	printf "[virusgatherer] transforming fastq to fasta'\n", $rid;
	for ( my $i=0; $i<=$#readFs; $i++ ){
		my $fqfile = $readFs[ $i ];  $fqfile =~ s/\.fastq/\.trim\.q$readQcut\.fastq/;
		my $fafile = $fqfile;        $fafile =~ s/\.fastq$/\.fasta/;
		my $cmd = sprintf "seqtk seq -a %s > %s", $fqfile, $fafile;
		`$cmd`;
	}
}


# -------- #
# clean up
# -------- #
if ( $tmpdir ne "" ){
	printf "[virusgatherer] moving results to original working directory %s\n", $orig_wdir;
	`cp     $tmpdir/$sraid/* $orig_wdir/`;
	`rm -rf $tmpdir/$sraid`;
}


# log end time and calculate runtime
$etime = time();
$rtime = $etime - $stime;
$ho = int( $rtime / 3600 );  $rtime -= $ho * 3600;
$mi = int( $rtime /   60 );  $rtime -= $mi *   60;
printf "[virusgatherer] done - runtime: %d h %d min %d sec\n", $ho, $mi, $rtime;

