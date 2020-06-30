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
my $resdir        = "/projects/p_sra/";
my $tmpdir        = sprintf "%s/lauber-2015-virusHunter", virusHGutils::get_workspace_path( "scratch" );
my $logtmpdir     = sprintf "%s/virushunter",             virusHGutils::get_workspace_path( "lustre"  );
#my $ncbirootASCP  = "130.14.29.30";
my $ncbirootFTP   = "130.14.250.7";
my $ebirootFTP    = "ftp.sra.ebi.ac.uk";
my $ebiIP         = "fasp.sra.ebi.ac.uk";
my $asperaroot    = "~/.aspera/connect";
my $ascpBandw     = 300;
my %sraIDsOld     = ();
my $sampleID      = "";
my $sraSize       = 0;
my $sraSizeLim    = 300; # 150? # [Gb]
my $delayDL       = 0; # [min]

# --------------- #
# usage and input
# --------------- #
if ( $#ARGV < 4 ){
die("
usage: virushunterHPC-1-data.pl <family> <project_id> <sra_run_id> <provider_ip> <provider> [options]\n
options:
\t-l=<float>\tdo not process runs larger than <l> Gb (default: $sraSizeLim)
\t-dly=<float>\tdelay consecutive downloads for <dly> minutes (default: $delayDL)
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
	if ( $option =~ /-l=(.*)/       ){ $sraSizeLim = $1; }
	if ( $option =~ /-dly=(\d+.*)/  ){ $delayDL    = $1; }
}


# -------------- #
# initialization
# -------------- #
init();


# -------------------------------- #
# download reads from SRA database
# -------------------------------- #
if ( $sraid =~ /^ERR*/ | $provider eq "ebi" ){
	dumpFromENA();
}else{
	dumpFromNCBI();
}


# --------- #
# functions
# --------- #
# init project (if started the first time)
sub init {
	$resdir    .= "/$family/results";
	$tmpdir    .= "/$family";
	$logtmpdir .= "/$family/$projID/log";
	# create directories and log file in case of a new project
	if ( -s "$resdir/$projID/log/completed.txt" > 0 ){
		open (OLD, "$resdir/$projID/log/completed.txt" );
		while ( my $line = <OLD> ){
			chomp($line);
			next if ( $line eq "" );
			$sraIDsOld{$line} = 1;
		}
		close(OLD);
		printf "[virushunter] project '%s' already existent and will be resumed\n", $projID;
	}
	else{
		open(EMPTY,">$resdir/$projID/log/completed.txt"); close(EMPTY);
		printf "[virushunter] new project '%s' initialized\n", $projID;
	}
	# was this SRA run already processed before?
	if ( exists $sraIDsOld{$sraid} ){
		printf "[virushunter] run '%s' already processed\n", $sraid;
		exit();
	}
} # end initProject

# download reads from NCBI
sub dumpFromNCBI {
	my ($stime, $etime, $rtime, $ho, $mi);
	# init new hits and runtime variables
	my %hits  = ();
	$stime = time();
	# download SRA run and transform to fasta
	printf "[virushunter] processing run '%s'; download from NCBI\n", $sraid;
	if ( ! -e "$tmpdir/$sraid.fasta" ){
		#my $dload  = "prefetch --ascp-path \"$asperaroot/bin/ascp|$asperaroot/etc/asperaweb_id_dsa.openssh\" --ascp-options \"-k 1 -QTr -l$ascpBandw"."m\" -v";
		my $dload  = "prefetch --max-size 100000000";
		   $dload .= " --output-file $tmpdir/$sraid.sra $sraid";
		printf "[virushunter] \t$sraid: starting download\n";
		printf "[virushunter] \t$dload\n";
		my $dload_out = `$dload`;
		printf STDERR "%s\n", $dload_out;
	}
	# failed?
	if ( ! -e "$tmpdir/$sraid.sra" or -s "$tmpdir/$sraid.sra" == 0 ){
		`rm -rf $tmpdir/$sraid.sra`;
		printf "[virushunter] \t$sraid: download unsuccessful - terminate\n";
		die "download unsuccessful - terminate\n";
	}
	printf "[virushunter] \t$sraid: data transfer completed\n";
	# get file size
	$sraSize = `ls -l $tmpdir/$sraid.sra | awk '{print \$5}'`;
	$sraSize /= ( 1024 * 1024 * 1024 );
	if ( $sraSize > $sraSizeLim ){
		printf "[virushunter] \t$sraid: run above the size limit of %.1f Gb - terminate\n", $sraSizeLim;
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
	open( DATA, ">$logtmpdir/$sraid-data.txt" );
	printf DATA "\n[virushunter-data]\n";
	printf DATA "Gb: %.4f\n", $sraSize;
	printf DATA "runtime: %d h %d min %d sec\n", $ho, $mi, $rtime;
	close(DATA);
	printf "[virushunter] \t$sraid: done - runtime: %d h %d min %d sec\n", $ho, $mi, $rtime;
	# log global statistics
	open(  STATS, ">>$tmpdir/virushunterHPC-$projID-data.txt" ) or die("$0: could not open file: $!\n");
	flock( STATS, 2 ) or die("$0: could not lock file: $!\n");
	printf STATS "%s\t%s\t%.4f\t%d\n", $sraid, $projID, $sraSize, $Tsec  if $Tsec > 0;
	close( STATS ) or die("$0: could not close file: $!\n");
	# sleep for a defined time period to delay next download
	sleep( $delayDL * 60 );
} # end dumpFromNCBI


# download reads from ENA
sub dumpFromENA{
	my ($stime, $etime, $rtime, $ho, $mi);
	# init new hits and runtime variables
	my %hits  = ();
	$stime = time();
	# download SRA run and transform to fasta
	printf "[virushunter] processing run '%s'; download from ENA\n", $sraid;
	if ( ! -e "$tmpdir/$sraid.fasta" ){
		my $dload = "$asperaroot/bin/ascp -i $asperaroot/etc/asperaweb_id_dsa.openssh -P33001 -k 1 -QTr -l$ascpBandw"."m ";
		my $srapath = "/".lc(substr($sraid,0,3))."/".substr($sraid,0,6)."/";
		if ( length($sraid) > 9 ){
			my $lastdig  = $sraid;
			   $lastdig  = chop($lastdig);
			   $srapath .= "00$lastdig/";
		}
		#$dload .= "era-fasp\@$providerIP:/vol1";
		$dload .= "era-fasp\@$ebiIP:/vol1";
		$dload .= $srapath.$sraid;
		$dload     .= " $tmpdir/$sraid.sra";
		printf "[virushunter] \t$sraid: starting download\n";
		#printf "$dload\n";
		my $ascpout = `$dload`;
		# try again ?
		if ( $ascpout =~ /ascp: failed to authenticate, exiting.*/  or  ! -e "$tmpdir/$sraid.sra" ){
			sleep(10);
			printf "[virushunter] \t$sraid: starting download, 2nd try\n";
			$ascpout = `$dload`;
		}
		# try using wget ?
		my $wgetout = "";
		if ( $ascpout =~ /ascp: failed to authenticate, exiting.*/  or  ! -e "$tmpdir/$sraid.sra" ){
			my $wd = cwd();  chdir( $tmpdir );
			printf "[virushunter] \t$sraid: starting download through FTP\n";
			$dload  = "wget --no-check-certificate -nv --tries=10 $ebirootFTP:/vol1";
			$dload .= $srapath.$sraid;			
			printf "$dload\n";
			$wgetout = `$dload`;
			if ( -e "$sraid" ){ `mv $sraid $sraid.sra`; }
			chdir( $wd );
		}
		# switch to the other provider and try a last wget ?
		if ( $wgetout =~ /.*ERROR 404: Not Found.*/  or  ! -e "$tmpdir/$sraid.sra" ){
			my $wd = cwd();  chdir( $tmpdir );
			printf "[virushunter] \t$sraid: switching provider and re-starting download through FTP\n";
			$dload  = "wget --no-check-certificate -nv --tries=10 $ncbirootFTP/sra/sra-instant/reads/ByRun/sra";
			$dload .= "/".substr($sraid,0,3)."/".substr($sraid,0,6)."/".$sraid."/$sraid.sra";
			printf "$dload\n";
			$wgetout = `$dload`;
			if ( -e "$sraid" ){ `mv $sraid $sraid.sra`; }
			chdir( $wd );
		}
		# failed?
		if ( ! -e "$tmpdir/$sraid.sra" or -s "$tmpdir/$sraid.sra" == 0 ){
			`rm -rf $tmpdir/$sraid.sra`;
			printf "[virushunter] \t$sraid: a last try via NCBI\n";
			dumpFromNCBI();
			return 0;
		}
		if ( ! -e "$tmpdir/$sraid.sra" or -s "$tmpdir/$sraid.sra" == 0 ){
			`rm -rf $tmpdir/$sraid.sra`;
			printf "[virushunter] \t$sraid: download unsuccessful - terminate\n";
			die "ascp and wget download unsuccessful - terminate\n";
		}
	}
	printf "[virushunter] \t$sraid: data transfer completed\n";
	# get file size
	$sraSize = `ls -l $tmpdir/$sraid.sra | awk '{print \$5}'`;
	$sraSize /= ( 1024 * 1024 * 1024 );
	if ( $sraSize > $sraSizeLim ){
		printf "[virushunter] \t$sraid: run above the size limit of %.1f Gb - terminate\n", $sraSizeLim;
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
	open( DATA, ">$logtmpdir/$sraid-data.txt" );
	printf DATA "\n[virushunter-data]\n";
	printf DATA "Gb: %.4f\n", $sraSize;
	printf DATA "runtime: %d h %d min %d sec\n", $ho, $mi, $rtime;
	close(DATA);
	printf "[virushunter] \t$sraid: done - runtime: %d h %d min %d sec\n", $ho, $mi, $rtime;
	# log global statistics
	open(  STATS, ">>$tmpdir/virushunterHPC-$projID-data.txt" ) or die("$0: could not open file: $!\n");
	flock( STATS, 2 ) or die("$0: could not lock file: $!\n");
	printf STATS "%s\t%s\t%.4f\t%d\n", $sraid, $projID, $sraSize, $Tsec  if $Tsec > 0;
	close( STATS ) or die("$0: could not close file: $!\n");
	# sleep for a defined time period to delay next download
	sleep( $delayDL * 60 );
} # end dumpFromENA


# download reads from SRA - !! deprecated !!
sub dumpFromSRA_ascp {
	my ($stime, $etime, $rtime, $ho, $mi);
	# init new hits and runtime variables
	my %hits  = ();
	$stime = time();
	# download SRA run and transform to fasta
	printf "[virushunter] processing run '%s'\n", $sraid;
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
			printf "[virushunter] \t$sraid: starting download, 2nd try\n";
			$ascpout = `$dload`;
		}
		# try using wget ?
		my $wgetout = "";
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
			printf "$dload\n";
			$wgetout = `$dload`;
			if ( -e "$sraid" ){ `mv $sraid $sraid.sra`; }
			chdir( $wd );
		}
		# switch to the other provider and try a last wget ?
		if ( $wgetout =~ /.*ERROR 404: Not Found.*/  or  ! -e "$tmpdir/$sraid.sra" ){
			my $wd = cwd();  chdir( $tmpdir );
			printf "[virushunter] \t$sraid: switching provider and re-starting download through FTP\n";
			if ( $provider eq "ncbi" ){
				$dload  = "wget --no-check-certificate -nv --tries=10 $ebirootFTP:/vol1";
				$srapath = "/".lc(substr($sraid,0,3))."/".substr($sraid,0,6)."/";
				if ( length($sraid) > 9 ){
					my $lastdig  = $sraid;
					   $lastdig  = chop($lastdig);
					   $srapath .= "00$lastdig/";
				}
				$dload .= $srapath.$sraid;			
			}else{
				$dload  = "wget --no-check-certificate -nv --tries=10 $ncbirootFTP/sra/sra-instant/reads/ByRun/sra";
				$dload .= "/".substr($sraid,0,3)."/".substr($sraid,0,6)."/".$sraid."/$sraid.sra";
			}
			printf "$dload\n";
			$wgetout = `$dload`;
			if ( -e "$sraid" ){ `mv $sraid $sraid.sra`; }
			chdir( $wd );
		}
		# failed?
		if ( ! -e "$tmpdir/$sraid.sra" or -s "$tmpdir/$sraid.sra" == 0 ){
			`rm -rf $tmpdir/$sraid.sra`;
			printf "[virushunter] \t$sraid: download unsuccessful - terminate\n";
			die "ascp and wget download unsuccessful - terminate\n";
		}
	}
	printf "[virushunter] \t$sraid: data transfer completed\n";
	# get file size
	$sraSize = `ls -l $tmpdir/$sraid.sra | awk '{print \$5}'`;
	$sraSize /= ( 1024 * 1024 * 1024 );
	if ( $sraSize > $sraSizeLim ){
		printf "[virushunter] \t$sraid: run above the size limit of %.1f Gb - terminate\n", $sraSizeLim;
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
	open( DATA, ">$logtmpdir/$sraid-data.txt" );
	printf DATA "\n[virushunter-data]\n";
	printf DATA "Gb: %.4f\n", $sraSize;
	printf DATA "runtime: %d h %d min %d sec\n", $ho, $mi, $rtime;
	close(DATA);
	printf "[virushunter] \t$sraid: done - runtime: %d h %d min %d sec\n", $ho, $mi, $rtime;
	# log global statistics
	open(  STATS, ">>$tmpdir/virushunterHPC-$projID-data.txt" ) or die("$0: could not open file: $!\n");
	flock( STATS, 2 ) or die("$0: could not lock file: $!\n");
	printf STATS "%s\t%s\t%.4f\t%d\n", $sraid, $projID, $sraSize, $Tsec  if $Tsec > 0;
	close( STATS ) or die("$0: could not close file: $!\n");
	# sleep for a defined time period to delay next download
	sleep( $delayDL * 60 );
} # end dumpAndMakeBlastDB

