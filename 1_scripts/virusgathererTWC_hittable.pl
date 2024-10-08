#!/usr/bin/perl

use warnings;
use strict;
use POSIX qw/floor/;
use Data::Dumper;


# input
if ( $#ARGV != 4 ){
die("
usage: virusgathererTWC-hittable.pl <parameters>\n
parameter:
\t<family>
\t<project>
\t<assembler>
\t<base_dir>
\t<is_SRA> [0|1]
\n");
}
my $family    = shift;
my $project   = shift;
my $assembler = shift;
my $basedir   = shift;
my $isSRA     = shift;
my $resdir    = "$basedir/$family/$project/results";
my $showDate  = 1;


# print header
if ( $isSRA ){
	printf "SRA_run\tSRA_sample\tSRA_study\thost_taxon\thost_taxid\tcontig_id\tcontig_len\tViralRefSeq_E\tViralRefSeq_ident\tViralRefSeq_aLen\tViralRefSeq_subject\tViralRefSeq_taxonomy\tdate_analyzed\n"; 
}else{
	printf "run_id\tcontig_id\tcontig_len\tViralRefSeq_E\tViralRefSeq_ident\tViralRefSeq_aLen\tViralRefSeq_subject\tViralRefSeq_taxonomy\tdate_analyzed\n"; 
}

# do for each data set
my @hitdirs = glob( "$resdir/*" );
foreach my $dir ( @hitdirs ){
	# results present ?
	my $resFile = "$dir/virusgatherer/genseedhmm-$assembler.fasta-vs-viral-blastx.tsv";
	next if ( ! -e "$resFile" );
	next if (   -s "$resFile" == 0 );
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
		$dateA = `stat --printf="%y\n" $resFile | awk '{print \$1}'`;
		chomp( $dateA );
	}
	
	# read virusgatherer results
	my %res = ();
	my ($txid, $txval, $accid, $blaststats);
	if ( -e "$resFile" ){
		open( RR, "<$resFile" ) or die( "Can't open final virusgatherer result file '$resFile': $!\n" );
		while ( my $line = <RR> ){
			chomp($line);  next if $line eq "";
			# parse blast result fields
			my @v       = split( /\t/, $line );
			$txid       = "";
			$txid       = $v[7]  if ( $#v > 6 );
			$accid      = $v[5];
			$blaststats = sprintf "%s\t%d\t%s\t%f\t%d\tacc:%s|%s", $v[0],$v[1],$v[2],$v[3],$v[4],$v[5],$v[6];
			# print for this hit
			if ( !$isSRA ){
				printf "%s\t%s",                 $sraid,                                   $blaststats;
			}else{
				printf "%s\t%s\t%s\t%s\t%s\t%s", $sraid, $srasam, $srastu, $taxon, $taxid, $blaststats;
			}
			# add taxonomy info
			$txval = getTaxForSeqId( $accid, "protein" );
			printf "\t%s", $txval;
			# finish line
			printf "\t%s\n", $dateA;
		}
	}
}

# function to get taxonomy for a nucleotide/protein accession/gi via eutils
sub getTaxForSeqId{
	# input
	my $sid  = shift;
	my $type = shift;
	# do the esearch
	my $esearch = "efetch -db $type -id $sid -format gpc";
	my $taxid   = "";
	my $taxa    = "";
	open( ESEARCH, "$esearch |" );
	while ( my $line = <ESEARCH> ){
		if ( $line =~ /<INSDSeq_taxonomy>(.*)<\/INSDSeq_taxonomy>/ ){
			my $tmp = $1;
			   $tmp =~ s/Viruses; //;
			my @tmp = split( /; /, $tmp );
			for ( my $i=$#tmp; $i>-1; $i--){  $taxa .= "|$tmp[$i]"; }
		}
		if ( $line =~ /<INSDQualifier_value>taxon:(\d+)<\/INSDQualifier_value>/ ){
			$taxid = $1;
		}
	}
	close(ESEARCH);
	# return result
	my $taxres = sprintf "taxid:%s%s", $taxid, $taxa;
	return( $taxres );
} # end getTaxForSeqId


