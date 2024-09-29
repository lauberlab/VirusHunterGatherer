#!/usr/bin/perl

# input
if ( $#ARGV != 0 ){
	die("\nusage: queryTax.pl <file_with_refseq_ids>\n\n");
}
my $idFile = shift;

# read RefSeq IDs from file
my @ids = ();
open( IDS, "<$idFile" ) or die( "Can't open file '$idFile': $!\n" );
while ( my $line = <IDS> ){
	chomp( $line );  next if $line eq "";
	$line =~ s/\|//g;
	push( @ids, $line );
}
close(IDS);

# get taxonomy for each RefSeq ID and send to STDOUT
foreach my $id ( @ids ){
	# use esearch to access taxonomy info
	my $esearch = "efetch -db nucleotide -id $id -format gpc";
	my $taxid   = "";
	my $taxa    = "";
	open( ESEARCH, "$esearch |" );
	while ( my $line = <ESEARCH> ){
		if ( $line =~ /<INSDSeq_accession-version>(.*)<\/INSDSeq_accession-version>/ ){
			my $returnid = $1;

			# Sanity check: Compare $returnid to $id
			# Why? Sometimes it returns not matching entry. Try "KNDV-Lp-2", for example
			# (should return empty result, but Bos taurus' sequence is returned O_o)
			if ( $returnid ne $id ){
				print STDERR "Error: Retrieved ID '$returnid' does not match the expected ID '$id'.\n";
				last;    # exit the loop
			}
		}
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
	# print to screen
	printf "%s\ttaxid:%s%s\n", $id, $taxid, $taxa;
}

