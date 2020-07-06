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
#use lib "/home/lauber/lib/perl/";
use warnings;
use strict;
use LWP::Simple qw/get/;
use POSIX qw/floor/;
use Cwd;
use Data::Dumper;

use lib "/home/lauber/lauber-2015-virusHunter";
use      virusHGutils;


# --------- #
# parameter
# --------- #
my $tmpdir        = sprintf "%s/lauber-2015-virusGatherer", virusHGutils::get_workspace_path( "scratch" );
#my $ncbirootASCP  = "130.14.29.30";
my $ncbirootFTP   = "130.14.250.7";
my $ebirootFTP    = "ftp.sra.ebi.ac.uk";
my $asperaroot    = "~/.aspera/connect";
my $ascpBandw     = 200;
my %sraIDsOld     = ();
my $sampleID      = "";
my $sraSize       = 0;
my $sraSizeLim    = 500; # 150? # [Gb]
my $useLustre     = 0;
my $userOutdir    = "";
my $userURL       = "";


# --------------- #
# usage and input
# --------------- #
if ( $#ARGV < 4 ){
die("
usage: virusgathererHPC-1-data.pl <family> <project_id> <sra_run_id> <provider_ip> <provider> [options]\n
options:
\t-l=<float>\tdo not process runs larger than <l> Gb (default: $sraSizeLim)
\t-lustre\t\tsave files on /lustre/ssd
\t-o=<string>\tuser-provided output directory
\t-url=<string>\tuser-provided URL to be used for download
\n");
}


# -------------------- #
# get input parameters
# -------------------- #
my $family     = shift;
my $projID     = shift;
my $sraid      = shift;
my $providerIP = shift;
my $provider   = shift;
while ( $#ARGV > -1 ){
	my $option = shift;
	if ( $option =~ /-l=(.*)/   ){ $sraSizeLim = $1; }
	if ( $option eq "-lustre"   ){ $useLustre  =  1; }
	if ( $option =~ /-o=(.*)/   ){ $userOutdir = $1; }
	if ( $option =~ /-url=(.*)/ ){ $userURL    = $1; }
}


# -------------- #
# initialization
# -------------- #
init();


# ----------------------------------------------- #
# download reads from SRA and blast against query
# ----------------------------------------------- #
if ( $userURL ne "" ){
	dumpData_URL();
}elsif( $sraid =~ /^ERR*/ ){
	dumpData_ENA();
}else{
	dumpData_NCBI();
}


# --------- #
# functions
# --------- #

# init project (if started the first time)
sub init {
	$tmpdir  = sprintf "%s/virusgatherer", virusHGutils::get_workspace_path( "lustre" )  if $useLustre;
	$tmpdir .= "/$family/$projID";
	$tmpdir  = $userOutdir  if $userOutdir ne "";
} # end initProject


# download reads from SRA
sub dumpData_NCBI {
	my ($stime, $etime, $rtime, $ho, $mi);
	# init new hits and runtime variables
	my %hits  = ();
	$stime = time();
	# download SRA run and transform to fasta
	printf "[virusgatherer] processing run '%s'; download via NCBI\n", $sraid;
	if ( ! -e "$tmpdir/$sraid.fasta" ){
		#my $dload  = "prefetch --ascp-path \"$asperaroot/bin/ascp|$asperaroot/etc/asperaweb_id_dsa.openssh\" --ascp-options \"-k 1 -QTr -l$ascpBandw"."m\" -v";
		my $dload  = "prefetch --max-size 100000000";
		   $dload .= " --output-file $tmpdir/$sraid.sra $sraid";
		printf "[virusgatherer] \t$sraid: starting download\n";
		printf "[virusgatherer] \t$dload\n";
		my $dload_out = `$dload`;
		printf STDERR "%s\n", $dload_out;
	}
	# failed?
	if ( ! -e "$tmpdir/$sraid.sra" or -s "$tmpdir/$sraid.sra" == 0 ){
		`rm -rf $tmpdir/$sraid.sra`;
		printf "[virusgatherer] \t$sraid: download unsuccessful - terminate\n";
		die "download unsuccessful - terminate\n";
	}
	printf "[virusgatherer] \t$sraid: data transfer completed\n";
	# get file size
	$sraSize = `ls -l $tmpdir/$sraid.sra | awk '{print \$5}'`;
	$sraSize /= ( 1024 * 1024 * 1024 );
	if ( $sraSize > $sraSizeLim ){
		printf "[virusgatherer] \t$sraid: run above the size limit of %.1f Gb - terminate\n", $sraSizeLim;
		`rm -rf $tmpdir/$sraid*`;
		die "SRA run too large - terminate\n";		
	}
	# runtime for this run
	$etime = time();
	$rtime = $etime - $stime;
	my $Tsec = $rtime;
	$ho = floor( $rtime / 3600 );  $rtime -= $ho * 3600;
	$mi = floor( $rtime /   60 );  $rtime -= $mi *   60;
	# log data runtime and SRA size
	open( DATA, ">$tmpdir/$sraid-data.txt" );
	printf DATA "\n[virusgatherer-data]\n";
	printf DATA "Gb: %.4f\n", $sraSize;
	printf DATA "runtime: %d h %d min %d sec\n", $ho, $mi, $rtime;
	close(DATA);
	printf "[virusgatherer] \t$sraid: done - runtime: %d h %d min %d sec\n", $ho, $mi, $rtime;
} # end dumpData_NCBI


# download reads from ENA
sub dumpData_ENA {
	my ($stime, $etime, $rtime, $ho, $mi);
	# init new hits and runtime variables
	my %hits  = ();
	$stime = time();
	# download SRA run and transform to fasta
	printf "[virusgatherer] processing run '%s'; download via ENA\n", $sraid;
	if ( ! -e "$tmpdir/$sraid.fasta" ){
		my $dload = "$asperaroot/bin/ascp -i $asperaroot/etc/asperaweb_id_dsa.openssh -P33001 -k 1 -QTr -l$ascpBandw"."m ";
		my $srapath = "/".lc( substr($sraid,0,3) )."/".substr($sraid,0,6)."/";
		if ( length($sraid) > 9 ){
			my $lastdig  = $sraid;
			$lastdig  = chop($lastdig);
			$srapath .= "00$lastdig/";
		}
		$dload .= "era-fasp\@$providerIP:/vol1";
		$dload .= $srapath.$sraid;		
		$dload .= " $tmpdir/$sraid.sra";
		printf "[virusgatherer] \t$sraid: starting download\n";
		#printf "$dload\n";
		my $ascpout = `$dload`;
		# try again ?
		if ( $ascpout =~ /ascp: failed to authenticate, exiting.*/  or  ! -e "$tmpdir/$sraid.sra" ){
			sleep(10);
			printf "[virusgatherer] \t$sraid: starting download, 2nd try\n";
			$ascpout = `$dload`;
		}
		# try using wget ?
		if ( $ascpout =~ /ascp: failed to authenticate, exiting.*/  or  ! -e "$tmpdir/$sraid.sra" ){
			my $wd = cwd();  chdir( $tmpdir );
			printf "[virusgatherer] \t$sraid: starting download through FTP\n";
			$dload  = "wget --no-check-certificate -nv --tries=10 $ebirootFTP:/vol1";
			$dload .= $srapath.$sraid;
			#printf "$dload\n";
			`$dload`;
			`mv $sraid $sraid.sra`;
			chdir( $wd );
		}
		if ( ! -e "$tmpdir/$sraid.sra" or -s "$tmpdir/$sraid.sra" == 0 ){
			`rm -rf $tmpdir/$sraid.sra`;
			printf "[virusgatherer] \t$sraid: a last try via NCBI\n";
			dumpData_NCBI();
			return 0;
		}
		if ( ! -e "$tmpdir/$sraid.sra" or -s "$tmpdir/$sraid.sra" == 0 ){
			`rm -rf $tmpdir/$sraid.sra`;
			printf "[virusgatherer] \t$sraid: download unsuccessful - terminate\n";
			die "ascp and wget download unsuccessful - terminate\n";
		}
	}
	printf "[virusgatherer] \t$sraid: data transfer completed\n";
	# get file size
	$sraSize = `ls -l $tmpdir/$sraid.sra | awk '{print \$5}'`;
	$sraSize /= ( 1024 * 1024 * 1024 );
	if ( $sraSize > $sraSizeLim ){
		printf "[virusgatherer] \t$sraid: run above the size limit of %.1f Gb - terminate\n", $sraSizeLim;
		`rm -rf $tmpdir/$sraid*`;
		die "SRA run too large - terminate\n";		
	}
	# runtime for this run
	$etime = time();
	$rtime = $etime - $stime;
	my $Tsec = $rtime;
	$ho = floor( $rtime / 3600 );  $rtime -= $ho * 3600;
	$mi = floor( $rtime /   60 );  $rtime -= $mi *   60;
	# log data runtime and SRA size
	open( DATA, ">$tmpdir/$sraid-data.txt" );
	printf DATA "\n[virusgatherer-data]\n";
	printf DATA "Gb: %.4f\n", $sraSize;
	printf DATA "runtime: %d h %d min %d sec\n", $ho, $mi, $rtime;
	close(DATA);
	printf "[virusgatherer] \t$sraid: done - runtime: %d h %d min %d sec\n", $ho, $mi, $rtime;
} # end dumpData_ENA


# download reads from SRA
sub dumpData_old {
	my ($stime, $etime, $rtime, $ho, $mi);
	# init new hits and runtime variables
	my %hits  = ();
	$stime = time();
	# download SRA run and transform to fasta
	printf "[virusgatherer] processing run '%s'\n", $sraid;
	if ( ! -e "$tmpdir/$sraid.fasta" ){
		my $dload = "$asperaroot/bin/ascp -i $asperaroot/etc/asperaweb_id_dsa.openssh -k 1 -QTr -l$ascpBandw"."m ";
		my $srapath = "";
		if ( $provider eq "ncbi" ){
			$dload .= "anonftp\@$providerIP:/sra/sra-instant/reads/ByRun/sra";
			$dload .= "/".substr($sraid,0,3)."/".substr($sraid,0,6)."/".$sraid."/$sraid.sra ";
		}else{
			$srapath = "/".lc(substr($sraid,0,3))."/".substr($sraid,0,6)."/";
			if ( length($sraid) > 9 ){
				my $lastdig  = $sraid;
				$lastdig  = chop($lastdig);
				$srapath .= "00$lastdig/";
			}
			$dload .= "era-fasp\@$providerIP:/vol1";
			$dload .= $srapath.$sraid;
		}
		$dload     .= " $tmpdir/$sraid.sra";
		printf "[virushunter] \t$sraid: starting download\n";
		#printf "$dload\n";
		my $ascpout = `$dload`;
		# try again ?
		if ( $ascpout =~ /ascp: failed to authenticate, exiting.*/  or  ! -e "$tmpdir/$sraid.sra" ){
			sleep(10);
			printf "[virusgatherer] \t$sraid: starting download, 2nd try\n";
			$ascpout = `$dload`;
		}
		# try using wget ?
		if ( $ascpout =~ /ascp: failed to authenticate, exiting.*/  or  ! -e "$tmpdir/$sraid.sra" ){
			my $wd = cwd();  chdir( $tmpdir );
			printf "[virushunter] \t$sraid: starting download through FTP\n";
			if ( $provider eq "ncbi" ){
				$dload  = "wget --no-check-certificate -nv --tries=10 $ncbirootFTP/sra/sra-instant/reads/ByRun/sra";
				$dload .= "/".substr($sraid,0,3)."/".substr($sraid,0,6)."/".$sraid."/$sraid.sra";
			}else{
				$dload  = "wget --no-check-certificate -nv --tries=10 $ebirootFTP:/vol1";
				$dload .= $srapath.$sraid;
			}
			#printf "$dload\n";
			`$dload`;
			`mv $sraid $sraid.sra`;
			chdir( $wd );
		}
		if ( ! -e "$tmpdir/$sraid.sra" or -s "$tmpdir/$sraid.sra" == 0 ){
			`rm -rf $tmpdir/$sraid.sra`;
			printf "[virusgatherer] \t$sraid: download unsuccessful - terminate\n";
			die "ascp and wget download unsuccessful - terminate\n";
		}
	}
	printf "[virusgatherer] \t$sraid: data transfer completed\n";
	# get file size
	$sraSize = `ls -l $tmpdir/$sraid.sra | awk '{print \$5}'`;
	$sraSize /= ( 1024 * 1024 * 1024 );
	if ( $sraSize > $sraSizeLim ){
		printf "[virusgatherer] \t$sraid: run above the size limit of %.1f Gb - terminate\n", $sraSizeLim;
		`rm -rf $tmpdir/$sraid*`;
		die "SRA run too large - terminate\n";		
	}
	# runtime for this run
	$etime = time();
	$rtime = $etime - $stime;
	my $Tsec = $rtime;
	$ho = floor( $rtime / 3600 );  $rtime -= $ho * 3600;
	$mi = floor( $rtime /   60 );  $rtime -= $mi *   60;
	# log data runtime and SRA size
	open( DATA, ">$tmpdir/$sraid-data.txt" );
	printf DATA "\n[virusgatherer-data]\n";
	printf DATA "Gb: %.4f\n", $sraSize;
	printf DATA "runtime: %d h %d min %d sec\n", $ho, $mi, $rtime;
	close(DATA);
	printf "[virusgatherer] \t$sraid: done - runtime: %d h %d min %d sec\n", $ho, $mi, $rtime;
} # end dumpData


# download read from SRA with given URL
sub dumpData_URL {
	my ($stime, $etime, $rtime, $ho, $mi);
	# init new hits and runtime variables
	my %hits  = ();
	$stime = time();
	# download SRA run and transform to fasta
	printf "[virusgatherer] processing run '%s'\n", $sraid;
	if ( ! -e "$tmpdir/$sraid.fasta" ){
		# using wget
		my $wd    = cwd();  chdir( $tmpdir );
		my $dload = "wget --no-check-certificate -nv --tries=10 $userURL";
		`$dload`;
		# rename downloaded file ?
		$userURL =~ /.*\/($sraid.*)$/;
		my $sraversion = $1;
		if ( $sraversion ne "$sraid.sra" ){
			`mv $sraversion $sraid.sra`;
		}
		chdir( $wd );
		# success ?
		if ( ! -e "$tmpdir/$sraid.sra" or -s "$tmpdir/$sraid.sra" == 0 ){
			`rm -rf $tmpdir/$sraid.sra`;
			printf "[virusgatherer] \t$sraid: download unsuccessful - terminate\n";
			die "wget download unsuccessful - terminate\n";
		}
	}
	printf "[virusgatherer] \t$sraid: data transfer completed\n";
	# get file size
	$sraSize = `ls -l $tmpdir/$sraid.sra | awk '{print \$5}'`;
	$sraSize /= ( 1024 * 1024 * 1024 );
	if ( $sraSize > $sraSizeLim ){
		printf "[virusgatherer] \t$sraid: run above the size limit of %.1f Gb - terminate\n", $sraSizeLim;
		`rm -rf $tmpdir/$sraid*`;
		die "SRA run too large - terminate\n";		
	}
	# runtime for this run
	$etime = time();
	$rtime = $etime - $stime;
	my $Tsec = $rtime;
	$ho = floor( $rtime / 3600 );  $rtime -= $ho * 3600;
	$mi = floor( $rtime /   60 );  $rtime -= $mi *   60;
	# log data runtime and SRA size
	open( DATA, ">$tmpdir/$sraid-data.txt" );
	printf DATA "\n[virusgatherer-data]\n";
	printf DATA "Gb: %.4f\n", $sraSize;
	printf DATA "runtime: %d h %d min %d sec\n", $ho, $mi, $rtime;
	close(DATA);
	printf "[virusgatherer] \t$sraid: done - runtime: %d h %d min %d sec\n", $ho, $mi, $rtime;
} # end dumpDataURL

