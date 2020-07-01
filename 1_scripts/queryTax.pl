#!/usr/bin/perl

# input
if ( $#ARGV != 0 ){
	die("\nusage: queryTax.pl <file_with_gis>\n\n");
}
my $giFile = shift;

# read GIs from file
my @gis = ();
open( GIS, "<$giFile" ) or die( "Can't open file '$giFile': $!\n" );
while ( my $line = <GIS> ){
	chomp( $line );  next if $line eq "";
	$line =~ s/\|//g;
	push( @gis, $line );
}
close(GIS);

# get taxonomy for each GI and send to STDOUT
foreach my $gi ( @gis ){
	# use esearch to access taxonomy info
	my $esearch = "efetch -db nucleotide -id $gi -format gpc";
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
	# print to screen
	printf "%s\ttaxid:%s%s\n", $gi, $taxid, $taxa;
}

