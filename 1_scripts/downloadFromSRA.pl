#!/usr/bin/perl

# input
if ( $#ARGV != 1 ){
	die( "\nusage: downloadFromSRA.pl <file_with_SRA_run_id> <output_directory>\n\n" );
}
my $srafile = shift;
my $outdir  = shift;

# verify output directory exists
`mkdir $outdir` if ( ! -d $outdir );

# get SRA ID from file
open( SID, "<$srafile" ) or die( "Can't read file '$srafile': $!\n" );
my $sraid = <SID>;  chomp($sraid);
close(SID);

# download
# when using -O, it creates additional $sraid direcory inside the output
#   direcory, where it puts the .sra. Therefore, let's specify the output file
#   instead (sra-tools=2.10.0):
`prefetch $sraid -o $outdir/$sraid.sra`;

# unpack
`fastq-dump -O $outdir -B --split-spot --skip-technical --readids --gzip --clip $outdir/$sraid.sra`;

# clean
`rm $outdir/$sraid.sra`;

