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
my $sra_file = "$outdir/$sraid.sra";

if (-e $sra_file) {
    my $validate_cmd = "vdb-validate $sra_file";
    my $validate_output = `$validate_cmd`;

    if ($? == 0) {
        print "Skipping prefetch. File validation successful.\n";
    } else {
        print "File validation failed. Executing prefetch...\n";
        `prefetch $sraid -o $sra_file`;
    }
} else {
    print "Executing prefetch...\n";
    `prefetch $sraid -o $sra_file`;
}

# unpack
`fastq-dump -O $outdir -B --split-spot --skip-technical --readids --gzip --clip $sra_file`;

# clean
`rm $sra_file`;

