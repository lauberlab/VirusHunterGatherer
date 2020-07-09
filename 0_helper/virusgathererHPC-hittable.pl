#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

use lib "/home/lauber/lauber-2015-virusHunter";
use      virusHGutils;

# parameters
my $ETH       = 1e-4;

# input
if ( $#ARGV != 2 ){
die("\nusage: virusgathererHPC-hittable.pl <family> <project> <assembler_used>\n\n");
}
my $family    = shift;
my $project   = shift;
my $assembler = shift;

# get project directories and read-in final contig info
my $tmpdir  = sprintf "%s/lauber-2015-virusGatherer/$family/$project", virusHGutils::get_workspace_path( "scratch" );
my %results = ();
foreach my $resdir ( glob( "$tmpdir/*" ) ){
	next if ( -f $resdir );
	$resdir =~ /.*\/([^\/]+)/;
	my $aid  =  $1;
	$resdir .= "/results-$assembler/final_result_dir";
	# get contig info
#	my $cfile   = "$resdir/final_contigs.fasta";
	my $cfile   = "$resdir/final_positive_contigs.fasta";
	my %contigs = ();
	if ( -e $cfile ){
		open( C, "<$cfile" ) or die("Can't read contig file '$cfile': $!\n");
		my $cid = "";
		while ( my $line = <C> ){
			chomp($line);  next if $line eq "";
			if ( $line =~ />(.*)/ ){
				$cid = $1;
				$contigs{$cid}  = "";
			}else{
				$contigs{$cid} .= $line;
			}
		}
		close(C);
	}
	$results{$aid} = ();
	$results{$aid}{'cN'} = scalar keys %contigs;
	$results{$aid}{'cL'} = 0;
	foreach my $cid ( keys %contigs ){
		$results{$aid}{'cL'} = length( $contigs{$cid} ) if ( length($contigs{$cid}) > $results{$aid}{'cL'} );
	}
	# get blast against virusDB info
	my $bfile = $cfile."-vs-virusDB-tblastx.tsv";
	$results{$aid}{'E'}     = "";
	$results{$aid}{'ident'} = "";
	$results{$aid}{'hitID'} = "";
	if ( -e $bfile ){
		open( B, "<$bfile" ) or die( "Can't read tblastx result file '$bfile': $!\n" );
		while ( my $line = <B> ){
			chomp($line);  next if $line eq "";
			my @v = split( /\t/, $line );
			next if $v[1] > $ETH;
			if ( $results{$aid}{'E'} eq "" ){      $results{$aid}{'E'}     .=   $v[1];
			}else{                                 $results{$aid}{'E'}     .= "/$v[1]"; }
			if ( $results{$aid}{'ident'} eq "" ){  $results{$aid}{'ident'} .=   $v[2];
			}else{                                 $results{$aid}{'ident'} .= "/$v[2]"; }
			if ( $results{$aid}{'hitID'} eq "" ){  $results{$aid}{'hitID'} .=   $v[3];
			}else{                                 $results{$aid}{'hitID'} .= "/$v[3]"; }
		}
		close(B);
	}
}

# print results to screen
printf "assembly_id\tcontig_num\tcontig_maxLen\ttblastx_E\ttblastx_ident\ttblastx_hit\n";
foreach my $aid ( sort { $results{$b}{'cL'} <=> $results{$a}{'cL'} } keys %results ){
	printf "%s\t%d\t%d\t%s\t%s\t%s\n", $aid, $results{$aid}{'cN'}, $results{$aid}{'cL'},
					   $results{$aid}{'E'}, $results{$aid}{'ident'}, $results{$aid}{'hitID'};
}

# summary
my ($cnt1, $cnt2 ) = (0, 0);
foreach my $aid ( keys %results ){
	if ( $results{$aid}{'cN'} > 0 ){  $cnt1++;
	}else{				  $cnt2++; }
}
printf STDERR "\n";
printf STDERR "finished:\t%d\n",     $cnt1;
printf STDERR "not finished:\t%d\n", $cnt2;
printf STDERR "\n";

