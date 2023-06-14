#!/usr/bin/perl

use warnings;
use strict;
use POSIX qw/floor/;
use Data::Dumper;


# input
if ( $#ARGV != 3 ){
die("
usage: virushunterTWC-hittable.pl <parameters>\n
parameter:
\t<family>
\t<project>
\t<base_dir>
\t<is_SRA> [0|1]
\n");
}
my $family   = shift;
my $project  = shift;
my $basedir  = shift;
my $isSRA    = shift;
my $Scut     = 1.0e-4;
my $RmPhage  = 1;
my $resdir   = "$basedir/$family/$project/results";
my $showDate = 1;


# print header
if ( $isSRA ){
	printf "SRA_run\tSRA_sample\tSRA_study\thost_taxon\thost_taxid\tnum_hits\tbest_E\tbest_query\tViralRefSeq_E\tViralRefSeq_ident\tViralRefSeq_aLen/sLen\tViralRefSeq_contigs\tViralRefSeq_subject\tViralRefSeq_taxonomy\tdate_analyzed\n"; 
}else{
	printf "run_id\tnum_hits\tbest_E\tbest_query\tViralRefSeq_E\tViralRefSeq_ident\tViralRefSeq_aLen/sLen\tViralRefSeq_contigs\tViralRefSeq_subject\tViralRefSeq_taxonomy\tdate_analyzed\n"; 
}

# do for each data set
my @hitdirs = glob( "$resdir/*" );
foreach my $dir ( @hitdirs ){
	next if ( ! -d "$dir" );
	next if ( ! -d "$dir/virushunter" );
	# SRA ID
	my $sraid = $dir;
	$sraid =~ s/.*\///;
	# get taxon from NCBI annotation
	my ($taxon, $taxid, $srastu, $srasam) = ("","","","");
	if ( $isSRA){
		my $command = "esearch -db sra -query $sraid | efetch -format native | xtract -pattern EXPERIMENT_PACKAGE -block SAMPLE -element SCIENTIFIC_NAME -element TAXON_ID -element \@accession -block STUDY -element \@accession";
		my $output = `$command`;
		# Remove leading and trailing whitespace from the output
		$output =~ s/^\s+|\s+$//g;
		# Split the output into tab-delimited fields
		my @fields = split(/\t/, $output);
		$taxon = $fields[0];
		$taxid = $fields[1];
		$srasam = $fields[2];
		$srastu = $fields[3];
	}

	# get date
	my $dateA = "";
	if ( $showDate == 1 ){
		$dateA = `stat --printf="%y\n" $dir/virushunter | awk '{print \$1}'`;
		chomp( $dateA );
	}
	
	# read virushunter results
	my %res = ();
	my $resFile = "$dir/virushunter/final.hits.tsv";
	if ( -e "$resFile" ){
		open( RR, "<$resFile" ) or die( "Can't open final virushunter result file '$resFile': $!\n" );
		while ( my $line = <RR> ){
			chomp($line);  next if $line eq "";
			# optionally, remove hits against phages
			if ( $RmPhage == 1 ){
				next if $line =~ /.*phage .*/i;
			}
			# print for this hit
			if ( !$isSRA ){
				printf "%s\t%s\t%s\n",                 $sraid,                                   $line, $dateA;
			}else{
				printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $sraid, $srasam, $srastu, $taxon, $taxid, $line, $dateA;
			}
		}
	}
}

