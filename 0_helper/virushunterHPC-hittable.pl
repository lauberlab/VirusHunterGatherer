#!/usr/bin/perl

use warnings;
use strict;
use POSIX qw/floor/;
use Data::Dumper;

use lib "/home/lauber/lauber-2015-virusHunter";
use      virusHGutils;


# input
if ( $#ARGV < 1 ){
die("
usage: virushunterHPC-hittable.pl <family> <project> [options]\n
options:
\t-s\t\tshow only a hit summary
\t-p\t\tdo not remove hits against phages
\t-ip\t\tanalysis is in progress (results still on /scratch)
\t-nosra\t\tsimply output directory IDs and hit tables
\t-date\t\tshow date of analysis for each SRA run
\t-d=yyyy-mm-dd\tonly consider hits found after specified date
\t-rsv=<string>\tuse reservation <rsv> instead of standard queue
\n");
}
my $family   = shift;
my $project  = shift;
my $Sonly    = 0;
my $RmPhage  = 1;
my $Scut     = 1.0e-4;
my $resdir   = "/projects/p_sra/$family/results";
my $tmpdir   = sprintf "%s/lauber-2015-virusHunter/$family",  virusHGutils::get_workspace_path( "scratch" );
my $inprgs   = 0;
my $deadline = "1981-09-15";
my $showDate = 0;
my $nosra    = 0;
my $reserv   = "";

# print only best E-value?
while ( $#ARGV > -1 ){
	my $option = shift;
	if ( $option eq "-s" ){       $Sonly    =  1; }
	if ( $option eq "-p" ){       $RmPhage  =  0; }
	if ( $option eq "-ip" ){      $inprgs   =  1; }
	if ( $option eq "-nosra"){    $nosra    =  1; }
	if ( $option =~ /-d=(.+)/ ){  $deadline = $1; }
	if ( $option eq "-date" ){    $showDate =  1; }
	if ( $option =~ /-rsv=(.*)/){ $reserv   = $1; }
}

# prepare for copy and extraction of archived data
printf STDERR "\nextracting archived data ...\n";
my $srun  = "srun --time=23:59:00 --quiet";
   $srun .= " --reservation=$reserv"  if( $reserv ne "" );
my $cmdir = `pwd`;  chomp( $cmdir );

# copy
mkdir( "$tmpdir/$project-hittable" );
chdir( "$tmpdir/$project-hittable" );
`$srun cp $resdir/$project/hits.tar.gz .`;

# extract
`$srun tar -zxvf hits.tar.gz`;
chdir( "$cmdir" );

# go through all hits of specified project
my $bestE   = 1000000;
my $bestID  = "";
my $ScutN   = 0;
my $StotN   = 0;
my $globdir = "$tmpdir/$project-hittable/hits";
   $globdir = "$tmpdir/$project"  if ( $inprgs == 1 );
my $hitdirs = `find $globdir -maxdepth 1 -newermt '$deadline'`;
my @hitdirs = split( /\n/, $hitdirs );
shift( @hitdirs )  if $hitdirs[0] =~ /.*\/hits$/;
shift( @hitdirs )  if $hitdirs[0] =~ /.*\/$project$/;

# print header
if ( $Sonly == 0 ){
	printf STDERR "\ncompiling hit table ...\n";
	if ( $nosra ){
		printf "ID\tnum_hits\tbest_E\tbest_query\tViralRefSeq_E\tViralRefSeq_ident\tViralRefSeq_aLen/sLen\tViralRefSeq_contigs\tViralRefSeq_subject\tViralRefSeq_taxonomy"; 
	}else{
		printf "SRA_run\tSRA_sample\tSRA_study\thost_taxon\thost_taxid\tnum_hits\tbest_E\tbest_query\tViralRefSeq_E\tViralRefSeq_ident\tViralRefSeq_aLen/sLen\tViralRefSeq_contigs\tViralRefSeq_subject\tViralRefSeq_taxonomy"; 
	}
}
if ( $showDate == 1 ){
	printf "\tdate_analyzed";
}
print "\n"  if ( $Sonly == 0 );

foreach my $dir ( @hitdirs ){
	next if ( ! -d $dir );
	# SRA ID
	my $sraid = $dir;
	$sraid =~ s/.*\///;
	# get taxon from NCBI annotation
	my ($taxon, $taxid, $srastu, $srasam) = ("","","","");
	if ( $Sonly == 0 and $nosra == 0 ){
		my $wget_cmd = "wget -q -O - 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=$sraid'";
		my $wget_res = "";
		my $wget_itr = 0;
		while ( $wget_res eq "" ){
			$wget_res = `$wget_cmd`;
			$wget_itr++;
			last if $wget_itr > 5;
			sleep( 12 ) if $wget_res eq "";
		}
		if ( $wget_res ne "" ){
			my @lin = split /\n/, $wget_res;
			my @hdr = split /,/,  $lin[0];
			my @val = split /,/,  $lin[1];
			for ( my $wget_i=0; $wget_i<=$#hdr; $wget_i++ ){
				$taxon  = $val[ $wget_i ] if ( $hdr[ $wget_i] eq "ScientificName" );
				$taxid  = $val[ $wget_i ] if ( $hdr[ $wget_i] eq "TaxID" );
				$srasam = $val[ $wget_i ] if ( $hdr[ $wget_i] eq "Sample" );
				$srastu = $val[ $wget_i ] if ( $hdr[ $wget_i] eq "SRAStudy" );
				$taxon =~ s/ /_/g;
			}
		}
	}
	my $dateA = "";
	if ( $showDate == 1 ){
		$dateA = `stat --printf="%y\n" $dir | awk '{print \$1}'`;
		chomp( $dateA );
	}
	
	# read virushunter results
	my %res = ();
	my $resFile = "$dir/final.hits.tsv";
	my $bE = 100000; my $bI="";
	if ( -e "$resFile" ){
		open( RR, "<$resFile" ) or die( "Can't open final virushunter result file '$resFile': $!\n" );
		while ( my $line = <RR> ){
			chomp($line);  next if $line eq "";
			my @v = split( /\t/, $line );
			if ( $v[1] < $bE ){
				$bE = $v[1];
				$bI = $v[2];
			}
			# fix virushunter bug not reporting taxonomy annotation
			#if ( $line !~ /.*taxid:.*/ ){
			#	$v[7] =~ /.*gi:(\d+)\|.*/;
			#	open(  GI, ">/scratch/lauber/tmp/queryTax-$$.txt" );
			#	printf GI "|%d|\n", $1;
			#	close( GI );
			#	my $taxraw = `/home/lauber/lauber-2015-virusHunter/virusutils-queryTax.pl /scratch/lauber/tmp/queryTax-$$.txt`;  chomp($taxraw);
			#	my @taxspl = split( /\t/, $taxraw );
			#	$line     .= sprintf "\t%s", $taxspl[1];
			#}
			# optionally, remove hits against phages
			if ( $RmPhage == 1 ){
				next if $line =~ /.*phage .*/i;
			}
			# print for this hit
			if ( $Sonly == 0 ){
				if ( $nosra ){
					printf "%s\t%s",                 $sraid,                                   $line;
				}else{
					printf "%s\t%s\t%s\t%s\t%s\t%s", $sraid, $srasam, $srastu, $taxon, $taxid, $line;
				}
			}
			# optionally, show date
			if ( $showDate == 1 ){
					printf "\t%s", 		     $dateA;
			}
			# end line
			printf "\n"  if ( $Sonly == 0 );
		}
	}
	# catalog some summary statistics
	if ( $Sonly == 1 ){
		if ( $bE < $bestE ){
			$bestE  = $bE;
			$bestID = $bI;
		}
		$ScutN++ if $bE < $Scut;
		$StotN++;
	}
}

# remove temporary files, if any
#`rm -rf /scratch/lauber/tmp/queryTax-$$.txt`;

# optionally print overall best E-value
if ( $Sonly == 1 ){
	printf STDERR "\ncompiling summary ...\n";
	my $sraN = `$srun tar -O -zxf $resdir/$project/log.tar.gz log/completed.txt | wc -l`;
	chomp($sraN);
	#$sraN =~ /^(\d+) .*/;
	#$sraN =  $1;
	# output result summary
	printf "best E-value: %.1e (%s)\n", $bestE, $bestID;
	if ( $ScutN > 0 ){  printf "hits with E-value < %.1e: %d\n", $Scut, $ScutN; }
	printf "total number of hits: %d",      $StotN;
	if ( $sraN  > 0 ){  printf " (%.1f%%)", $StotN / $sraN * 100; } 
	print  "\n";
}

# clean up
print STDERR "\ncleaning up ...\n";
`$srun rm -rf $tmpdir/$project-hittable`;

# END
printf STDERR "\ndone\n\n";

