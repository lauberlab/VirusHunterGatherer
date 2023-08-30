#!/usr/bin/perl

# Genseed-HMM v. 1.0.6, 2015-11-10
# - first release version

use warnings;
use strict;
use Storable;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use Carp ();

# in case of problems, uncomment the following two lines for very detailed warning/error messages
#local $SIG{__WARN__} = \&Carp::cluck;
#local $SIG{__DIE__} = \&Carp::confess;

# turn on buffer autoflush
$| = 1;

# Capture the commands given to Genseed-HMM in command line for printing it after in genseed.log
my $script_command = $0;
foreach (@ARGV) {
    $script_command .= /\s/ ?   " \"" . $_ . "\""
                    :           " "   . $_;
}
my $version = "1.0.6";
my $last_update = "2015-11-10";
my $help_print = "GenSeed-HMM version $version ($last_update)

Usage: genseed-hmm.pl -seed <seed filename> -db <database filename> -assembler <assembler name> -ext_seed_size <integer>  <optional parameters> 

       or

       genseed-hmm.pl -conf <configuration file>

Mandatory parameters:
-seed                         Seed file name (containing profile HMM, DNA or protein sequence)
-db                           Read database file name 
-assembler                    Select your assembler (ABySS|Velvet|SOAPdenovo|Newbler|CAP3)  
-ext_seed_size <integer>      Length of extension seeds (contig ends to be used as seeds in the next round)

OPTIONAL PARAMETERS:
-conf                         Configuration file
-duplicate_headers <yes|no>   The database contains duplicate headers, e.g. due to paired reads (default = no)
-exp_direction                Expansion direction 
                                 left = expand only from the 5' end of the contig(s)
                                 right = expand only from the 3'end of the contig(s)
                                 both = expand to both ends of the contig(s) (default)
-threads <integer>            Number of processors to be used in GenSeed-HMM and third-party programs (default = 1)
-kmer <integer>               Kmer length to be used in Velvet, ABySS, and SOAPdenovo assemblers (default = 31)
-max_contig_length <integer>  Maximum contig length (default = no limit)
-max_number_rounds <integer>  Maximum number of rounds (default = no limit)
-min_hsp_length <integer>     Minimum length of the BLAST's HSP in percentage, relative of the seed sequence length (default = 0)
-output                       Output directory (default = genseed-hmm_dir#)
-clean <yes|no>               Remove all intermediate files, generate only final assembly files (default = yes)                              
-mapping <yes|no>             Map all selected reads onto positive contigs in the last round and generate a SAM file (default = yes)
-quality_inspector <integers> Trim contig ends to remove low quality bases (default = no trimming)
                                 first integer  = Length of the contig quality verification window
                                 second integer = Minimum quality value for the sliding window
                                 third integer  = Maximum percentage of bad quality values in the sliding window
-offset                       Offset for fastq conversion to fasta + qualities. Use 33 for Sanger (default) or 64 for Solexa
-align_threshold              Minimum ratio, in percentage, of alignment and previous contig lengths (default: 95)
-use_qual <yes|no>            Create quality files (if database of reads contains qualities) for use by the assembler (default = no)
-help                         Print this help and exit

IMPORTANT: the following parameters MUST be defined within double quotes in the command line
-blastn_parameters            blastn parameters (example: \"-dust no -perc_identity 90\")
-tblastn_parameters           tblastn parameters (example: \"-seg no -evalue 1\")
-hmmsearch_parameters         hmmsearch parameters (example: \"-T 1 --cpu 6\")
-assembler_parameters         assembler parameters\n";

#
# mandatory fields
#
my $seed;
my $db;
my $assembler;
my $assembler_name;
my $next_seed_size;

#
# parameters related to the database
#
my $duplicate_headers = 'no';
my $db_info->{'db_type'} = '';
my $fasta_qual = ' ';

#
# number of threads
#
my $threads = 1;

#
# parameters related to the aligners
#
my $blastn_parameters = '';
my $tblastn_parameters = '';
my $hmmsearch_parameters = '';

#
# parameters related to the assemblers
#
my $kmer = 31;
my $assembler_parameters = "";
my @quality_inspector;

#
# parameters exclusively related to GenSeed-HMM
#
my $max_number_rounds;
my $min_hsp_length = 0;
my $max_contig_length = 1000000;
my $output = 'genseed-hmm_dir';
my $clean = 'yes';
my $mapping = 'yes';
my $exp_direction = 'both';
my $help;
my $offset = 33;
my $align_threshold = 95;
my $redundancy = 10;
my $use_qual = 0;

#
# parameters related to the configuration file
#
my @arguments;
my %config;
my $missingArgument=0;
my $conf;

my $optret = GetOptions ("conf=s"                  => \$conf,
			 "seed=s"                  => \$seed,
			 "db=s"                    => \$db,
			 "duplicate_headers=s"     => \$duplicate_headers,
			 "threads=i"               => \$threads,
			 "exp_direction=s"         => \$exp_direction,
			 "ext_seed_size=i"         => \$next_seed_size,
			 "max_number_rounds=i"     => \$max_number_rounds,
			 "min_hsp_length=i"        => \$min_hsp_length,
			 "max_contig_length=i"     => \$max_contig_length,
			 "align_threshold=i"       => \$align_threshold,
			 "output=s"                => \$output,
			 "clean=s"                 => \$clean,
			 "blastn_parameters=s"     => \$blastn_parameters, 
			 "tblastn_parameters=s"    => \$tblastn_parameters,
			 "hmmsearch_parameters=s"  => \$hmmsearch_parameters, 
			 "assembler=s"             => \$assembler, 
			 "assembler_parameters=s"  => \$assembler_parameters,
			 "kmer=i"                  => \$kmer, 
			 "redundancy=i"            => \$redundancy, 
			 "mapping=s"               => \$mapping, 
			 "quality_inspector=i{3}"  => \@quality_inspector,
			 "use_qual=s"              => \$use_qual,
			 "offset=i"                => \$offset,
			 "help"                    => \$help,
    );

unless($optret) {
    die "\nERROR: problem found with options, please check above.\n";
}

if ($help) { 
    die $help_print; 
}

#
# variables for ID/offset tables created by indexing
#
my $id_table = {};
my $qid_table = {};

#
# if configuration file is not defined, check if the mandatory arguments are defined (seed file and database file)
#
if (not defined $conf) {
    if ( (!($seed)) or (!-e($seed)) ) { 
	print "ERROR: Missing mandatory argument -seed, or file not found.\n\n";
	die $help_print;
    }
    elsif ( (!($db)) or (!-e($db)) ) {
	print "ERROR: Missing mandatory argument -db, or file not found.\n\n";
	die $help_print;
    }
    elsif(!$assembler) {
	print "ERROR: Missing mandatory argument -assembler.\n\n";
	die $help_print;
    }
    elsif(!$next_seed_size) {
	print "ERROR: Missing mandatory argument -ext_seed_size.\n\n";
	die $help_print;
    }
}

if ($conf) {
    print STDERR "Configuration file specified, command line options will be overridden.\n";

    open(CONFIG, "< $conf") or die("ERROR: Problem opening configuration file $conf: $!\n");

    my $configLine;
    while ($configLine = <CONFIG>) {
        $configLine =~ s/^\s+//;
        $configLine =~ s/\s+\Z//;
        $configLine =~ s/\s+\=/\=/;
        $configLine =~ s/\=\s+/\=/;
        if ($configLine =~ /^\#/ || !($configLine =~ /(.)\=(.)/)) {
            next;
        }
        chomp $configLine;

	if ($configLine =~ m/(.+?)=(.+)/) {
	    $config{$1} = $2;
	}
    }
    close(CONFIG);

    #
    # Mandatory arguments
    #
    if (!($config{"seed"}) or (!-e($config{"seed"}))) {
	$missingArgument = 1;
	print "Missing mandatory configuration argument or file not found: seed.\n";
    }
    else {
	$seed = $config{"seed"};
    }
    if (!($config{"db"}) or (!-e($config{"db"}))) {
	$missingArgument = 1;
	print "Missing mandatory configuration argument or file not found: db.\n";
    }
    else {
	$db = $config{"db"};
    }
    if(!$config{"assembler"}) {
	$missingArgument = 1;
        print "Missing mandatory configuration argument: assembler.\n";
    }
    else {
	$assembler = $config{"assembler"};
    }
    if(!$config{"ext_seed_size"}) {
	$missingArgument = 1;
	print "Missing mandatory configuration argument: ext_seed_size.\n";
    }
    else {
	$next_seed_size = $config{"ext_seed_size"};
    }
    
    if ($missingArgument) {
	die "\nERROR: Cannot run genseed-hmm.pl, mandatory configuration argument(s) missing (see above).\n";
    }

    #
    # Optional parameters
    #
    if (defined($config{"duplicate_headers"})) {
	$duplicate_headers = $config{"duplicate_headers"};
    }
    if (defined($config{"threads"})) {
	$threads = $config{"threads"};
    }
    if (defined($config{"kmer"})) {
	$kmer = $config{"kmer"};
    }
    if (defined($config{"max_number_rounds"})) {
	$max_number_rounds = $config{"max_number_rounds"};
    }
    if (defined($config{"min_hsp_length"})) {
	$min_hsp_length = $config{"min_hsp_length"};
    }
    if (defined($config{"max_contig_length"})) {
	$max_contig_length = $config{"max_contig_length"};
    }
    if (defined($config{"blastn_parameters"})) {
	$blastn_parameters = $config{"blastn_parameters"};
    }
    if (defined($config{"tblastn_parameters"})) {
	$tblastn_parameters = $config{"tblastn_parameters"};
    }
    if (defined($config{"hmmsearch_parameters"})) {
	$hmmsearch_parameters = $config{"hmmsearch_parameters"};
    }
    if (defined($config{"assembler_parameters"})) {
	$assembler_parameters = $config{"assembler_parameters"};
    }
    if (defined($config{"align_threshold"})) {
	$align_threshold = $config{"align_threshold"};
    }
    if (defined($config{"redundancy"})) {
	$redundancy = $config{"redundancy"};
    }
    if (defined($config{"exp_direction"})) {
	$exp_direction = $config{"exp_direction"};
    }
    if (defined($config{"mapping"})) {
	$mapping = $config{"mapping"};
    }
    if (defined($config{"clean"})) {
	$clean = $config{"clean"};
    }
    if (defined($config{"use_qual"})) {
	$use_qual = $config{"use_qual"};
    }
    if (defined($config{"quality_inspector"})) {
	@quality_inspector = split(/\s/,$config{"quality_inspector"});
    }
    if (defined($config{"output"})) {
	$output = output_dir_name($config{"output"});
    }
    else {
	$output = output_dir_name($output);
    }
}
else {
    $output = output_dir_name($output);
}

$assembler_name = $assembler;

my $start_time = time;

#
# add multithreading information to search program parameters if it was not done by the user
#
$blastn_parameters .= " -num_threads $threads" unless $blastn_parameters =~ m/num_threads/;
$tblastn_parameters .= " -num_threads $threads" unless $tblastn_parameters =~ m/num_threads/;
$hmmsearch_parameters .= " --cpu $threads" unless $hmmsearch_parameters =~ m/cpu/i;

#
# store seed information
#
my $seed_information;
$seed_information->{'name'} = $seed;
$seed_information->{'type'} = seed_type_verify($seed_information->{'name'});
$seed_information->{'original_name'} = $seed_information->{'name'};
$seed_information->{'original_type'} = $seed_information->{'type'};

#
# Check programs
#
my @program_list = qw[hmmsearch tblastn blastn splitter makeblastdb bowtie2-build bowtie2 transeq sfffile sffinfo];
my @assembler_list = qw[abyss SOAPdenovo-31mer SOAPdenovo-63mer SOAPdenovo-127mer velveth velvetg runAssembly cap3];
my %error_hash;

print "Checking programs... \n";

foreach my $program_name (@program_list) {
    my $check = `which $program_name 2> /dev/null`;

    if (!$check) {
        print STDERR "\tProgram $program_name not found in path.\n";
        $error_hash{$program_name}++;
    }
}
foreach my $program_name (@assembler_list) {
    my $check = `which $program_name 2> /dev/null`;

    if (!$check && $program_name =~ /$assembler_name/i) {
        print STDERR "\tERROR: Program $program_name (assembler) not found in path but selected by user ($assembler_name).\n";
        $error_hash{$program_name}++;
    }
}
print "Done.\n";

#
# BLAST+ programs
#
if ($error_hash{'blastn'} || $error_hash{'makeblastdb'}) {
    print "ERROR: BLAST+ package programs not found. Cannot run GenSeed-HMM without BLAST+\n";
    exit;
}

#
# Check assemblers
#
if(($error_hash{'runAssembly'} || $error_hash{'splitter'}) && $assembler_name =~ m/Newbler/i) {
    print "ERROR: Could not find splitter or Newbler in the command path. Please correct or select another assembler.\n";
    exit;
}
if($error_hash{'cap3'} && $assembler_name =~ m/CAP3/i) {
    print "ERROR: Could not find cap3 program in the command path. Please correct or select another assembler.\n";
    exit;
}
if(($error_hash{'velveth'} || $error_hash{'velvetg'}) && $assembler_name =~ m/Velvet/i) {
    print "ERROR: Could not find velveth or velvetg in the command path. Please correct or select another assembler.\n";
    exit;
}
if ($error_hash{'abyss'} && $assembler_name =~ m/ABySS/i) {
    print "ERROR: Could not find abyss in the command path. Please correct or select another assembler.\n";
    exit;
} 
if (($error_hash{'SOAPdenovo-31mer'} || $error_hash{'SOAPdenovo-63mer'} || $error_hash{'SOAPdenovo-127mer'}) && $assembler_name =~ m/SOAPdenovo/i) {
    print "ERROR: Could not find SOAPdenovo in the command path. Please correct or select another assembler.\n";
    exit;
}

#
# Check other programs
#
if(($error_hash{'sffinfo'} || $error_hash{'sfffile'}) && $assembler_name =~ m/Newbler/i && $db =~ /\.sff$/i) {
    print "ERROR: Could not find sffinfo or sfffile, needed for Newbler assembly of SFF files, in the command path. Please correct, or select another assembler.\n";
    exit;
}

#
# Can GenSeed-HMM use protein seeds?
#
if ($error_hash{'tblastn'} && $seed_information->{'original_type'} =~ /Protein/) {
    print "ERROR: GenSeed-HMM can not start sequence reconstruction using a protein seed without the tblastn program.\n";
    exit;
}

#
# Can GenSeed-HMM use Profile HMM seeds?
#
if($error_hash{'hmmsearch'} && $seed_information->{'original_type'} =~ /Profile HMM/) {
    print "ERROR: GenSeed-HMM can not start the sequence reconstruction using a profile HMM seed without the hmmsearch program.\n";
    exit;
}

#
# Can GenSeed-HMM perform the final mapping?
#
if (($error_hash{'bowtie2-build'} || $error_hash{'bowtie2'}) && $mapping =~ m/yes/i) {
    print "WARNING: GenSeed-HMM can not perform the final mapping without the bowtie2 program. Turning off the -mapping parameter.\n";
    $mapping = 'no';
}

#
# create output directory
#
make_output_dir();

#
# creating genseed-hmm.log
#
my $log = "$output/genseed-hmm.log";
open(my $log_file_handle, ">$log") or die "ERROR: Could not create log file $log (548): $!\n";

#
# define file extension, contingent on read database type
#
my $ext;

# Check if the given db is fasta, fastq or SFF
#
print "Checking read database type... ";
print_current_time("Checking read database type... ");

unless( open(DATA, "<$db") ) {
    print STDERR "ERROR: Could not open database file $db (561): $!\n\n";
    exit;
}

while (my $line = <DATA>) {
    # check if DB is FASTA
    if ($line =~ /^>/) {
	$ext = "fasta";
	$db_info->{'db_type'} = "fasta"; 
	$fasta_qual = $db . ".qual" if (-e "$db\.qual");
	last;
    }
    # check if DB is FASTQ
    elsif($line =~ /^@/) {
	$ext = "fastq";
	$db_info->{'db_type'} = "fastq";            
	last;
    }
    else {
        $db_info->{'db_type'} = "sff";
	$ext = "sff";
	last;
    }
}

close DATA;

print "Detected: $db_info->{'db_type'}\n";
print_current_time("Detected: $db_info->{'db_type'}\n");

#
# print log information in genseed-hmm.log
#
my $date = `date`;
print $log_file_handle "# GenSeed-HMM :: a seed-driven progressive assembly program\n";
print $log_file_handle "# GenSeed-HMM is an Open Source program. Contributions for improving the program are welcome.\n";
print $log_file_handle "# Freely distributed under the GNU General Public License v. 3\n\n";
print $log_file_handle "#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
print $log_file_handle "Program called as:         $script_command\n";
print $log_file_handle "#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
print $log_file_handle "Program version:           $version\n";
print $log_file_handle "Date:                      $date";
print $log_file_handle "Seed:                      $seed_information->{original_name}\n";
print $log_file_handle "Seed type:                 $seed_information->{original_type}\n";
print $log_file_handle "Database:                  $db\n";
print $log_file_handle "Database type:             $db_info->{'db_type'}\n";
print $log_file_handle "FASTA quality file:        $fasta_qual\n" if ($db_info->{'db_type'} =~ m/fasta/i) && (-e "$db\.qual");
print $log_file_handle "Duplicate headers name:    $duplicate_headers\n";
print $log_file_handle "Expansion direction:       $exp_direction\n";
print $log_file_handle "Number of threads:         $threads\n";
print $log_file_handle "Assembler used:            $assembler_name\n";
print $log_file_handle "Kmer length:               $kmer\n" if ($assembler_name =~ m/(abyss)|(velvet)|(soapdenovo)/i);
print $log_file_handle "Assembler parameters:      $assembler_parameters\n";
print $log_file_handle "Extension seed size:       $next_seed_size\n";
print $log_file_handle "Output directory:          $output\n";
print $log_file_handle "blastn parameters:         $blastn_parameters\n";
print $log_file_handle "tblastn parameters:        $tblastn_parameters\n" if ($seed_information->{'original_type'} =~ m/Protein/);
print $log_file_handle "hmmsearch parameters:      $hmmsearch_parameters\n" if ($hmmsearch_parameters);
print $log_file_handle "Maximum contig length:     $max_contig_length\n" if ($max_contig_length != 1000000);
print $log_file_handle "Minimum HSP length:        $min_hsp_length\n";
print $log_file_handle "Clean output:              $clean\n";
print $log_file_handle "Maximum number of rounds:  $max_number_rounds\n" if (defined $max_number_rounds);
print $log_file_handle "Mapping:                   $mapping\n";
if (@quality_inspector) { 
    print $log_file_handle "Quality window:            @quality_inspector\n" if ($assembler_name =~ m/(cap3)|(newbler)/i)
}
else {
    print $log_file_handle "Quality window:            no quality trimming\n" if ($assembler_name =~ m/(cap3)|(newbler)/i)
}
print $log_file_handle "#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n";

if ($seed_information->{'original_type'} =~ m/DNA/) {
    print "\nSeed type: DNA\n";
}
elsif ($seed_information->{'original_type'} =~ m/Protein/) {
    print "\nSeed type: Protein\n";
}
elsif ($seed_information->{'original_type'} =~ m/Profile HMM/) {
    print "\nSeed type: Profile HMM\n";
}
print "Assembler: $assembler_name\n\n";


# If the DB is FASTA, continue. If it's FASTQ, check if is already formatted, and if it's not, convert it to FASTA
if($db_info->{'db_type'} =~ m/fasta/i) {
    $ext = "fasta"; 
} 
elsif($db_info->{'db_type'} =~ m/fastq/i) {
    if(-e $db.".fasta" && -e $db.".fasta.qual") {
	print "\nDatabase has already been converted from FASTQ to FASTA. Setting -db_type to fasta.\n";
	$ext = "fasta";
	$db_info->{'db_type'} = "fasta";
	$db = $db.".fasta";
	$fasta_qual = $db.".qual";
    }
    else {
	#
	# convert FASTQ to FASTA and continue as if input database had been FASTA all along
	#
	$db_info->{'db_type'} = "fasta";
	$db = fastq2fastaseq($db);
	$ext = "fasta";
	$fasta_qual = "$db\.qual";
    }
}
elsif($db_info->{'db_type'} =~ m/sff/i) {
    $ext = "sff";
}

#
# format headers if necessary and create blast database
#
if ($duplicate_headers =~ m/yes/i) {
    $db = format_multifasta_headers($db);

    if($fasta_qual ne ' ') {
	$fasta_qual = format_multifasta_headers($fasta_qual, "$db\.qual");
    }
}

#
# SFF file check
#
if($db =~ /\.sff$/i && $assembler_name =~ m/newbler/i) { 
    if($db_info->{'db_type'} !~ m/sff/i) { 
        print "WARNING: Option db_type was set to $db_info->{'db_type'} but assembler_name is $assembler_name and sff was given ($db);\nSwitching -db_type to sff.\n\n";
	$db_info->{'db_type'} = 'sff'; 
    }
}

#
# run makeblastdb to format BLAST database
#
run_makeblastdb($db);

$db_info->{'name'} = $db;
if ($db_info->{'db_type'} eq 'sff') { $db_info->{'name'} = "SFFextracted.fasta"; }

#
# if seed type is HMM, create translated database with transeq (EMBOSS)
#
if($seed_information->{'type'} =~ m/Profile HMM/) {
    my $name;
    my $extension = "";

    if ($db_info->{'name'} =~ m/(.*)\.([\w]+)$/){
	$name = $1;
	$extension = $2;
    }
    else {
	$db_info->{'name'} = $db;
    }

    my $transeq_db = $name."_transeq_result.".$extension;

    if (!-e $transeq_db) {
	run_transeq($db_info->{'name'}, $transeq_db, $threads);
    }
    $db_info->{'transeq'} = $transeq_db;
}

#
# start progressive assembly process
#

print $log_file_handle "\nStarting progressive assembly process... \n\n";
print "\nStarting progressive assembly process... \n\n";

my ($all_recruited_reads, $final_contigs) = progressive_assembly_routine();

print $log_file_handle "Finishing progressive assembly process...\n\n";
print $log_file_handle "###     Final check     ###\n";
print "Finishing progressive assembly process...\n\n";
print "###     Final check     ###\n";
print "Recruited ", scalar keys %$all_recruited_reads, " reads in the progressive assembly.\n";
print $log_file_handle "Recruited ", scalar keys %$all_recruited_reads, " reads in the progressive assembly.\n";

#
# remove all intermediate files if clean = yes 
#
if ($clean =~ m/yes/i) {
    unlink "$output/assembler.log" if (-e "$output/assembler.log");
    unlink "$output/check_progress.txt" if (-e "$output/check_progress.txt");
    unlink "$output/input_sequences.$ext" if (-e "$output/input_sequences.$ext");
    unlink "$output/input_sequences.$ext\.qual" if (-e "$output/input_sequences.$ext\.qual");
    unlink "$output/seed_file.fasta" if (-e "$output/seed_file.fasta");
    unlink "$output/tabular.results" if (-e "$output/tabular.results");
    unlink "$output/tabular_temp.results" if (-e "$output/tabular_temp.results");
    unlink "SFFextracted.fasta" if (-e "SFFextracted.fasta");
    unlink "$db\.fasta" if (-e "$db\.fasta");
    unlink "$db\.ref" if (-e "$db\.ref");
    unlink "$fasta_qual\.ref" if (-e "$fasta_qual\.ref");
}

#
# verify if the last contig(s) have the original seed
#
if (%$final_contigs) {
    open(FINAL_CONTIG, ">>$output/final_result_dir/final_contigs.fasta") or
	die "ERROR: Could not write to file $output/final_result_dir/final_contigs.fasta (669): $!\n";
    foreach my $name (keys %$final_contigs) {
	print FINAL_CONTIG ">$name\n";
	print FINAL_CONTIG "$final_contigs->{$name}\n";
    }
    close FINAL_CONTIG;
}

print "Checking if seed is present in the final contigs... \n";

my $final_positive_contigs = {};

check_progress($seed_information->{'original_name'}, "$output/final_result_dir/final_contigs.fasta", $seed_information->{'original_type'}, $output, $final_positive_contigs, "last");

if (!%$final_positive_contigs) {
    print "No contig is positive for the original seed file. ABORTING!\n";
    print_current_time("No contig is positive to the original seed file. ABORTING!");
    system "mv $output/check_progress.txt $output/check_progress_dir/final_check.txt" if (-e "$output/check_progress.txt" && $clean =~ m/no/i);
    unlink "$output/454AllContigs.qual" if (-e "$output/454AllContigs.qual");
    unlink "$output/contigs_next_round.fasta" if (-e "$output/contigs_next_round.fasta");
    unlink "$output/makeblastdb.log" if (-e "$output/makeblastdb.log" && $clean !~ m/no/i);
    exit;
}

#
# print the final positive configs file
#
print "Done.\nReconstructed ", scalar keys %$final_positive_contigs, " positive contig(s). Check $output/final_result_dir/final_positive_contigs.fasta\n";
print $log_file_handle "Reconstructed ", scalar keys %$final_positive_contigs, " positive contig(s). Check $output/final_result_dir/final_positive_contigs.fasta\n";
print "----------------------------------------------\n";
print $log_file_handle "----------------------------------------------\n";
print "Length of final positive contig(s):\n";
print $log_file_handle "Length of final positive contig(s):\n";

open(FINAL, ">$output/final_result_dir/final_positive_contigs.fasta") or
    die "ERROR: Could not create file $output/final_result_dir/final_positive_contigs.fasta (705): $!\n";

my $count  = 1;

foreach my $name (sort { length $final_positive_contigs->{$b} <=> length $final_positive_contigs->{$a} } keys %$final_positive_contigs) {
    print FINAL ">Contig-$count\n";
    print "-Contig-$count with ", length $final_positive_contigs->{$name}, " bp\n";
    print $log_file_handle "-Contig-$count with ", length $final_positive_contigs->{$name}, " bp\n";
    $count++;

    for (my $position = 0; $position < length $final_positive_contigs->{$name}; $position += 60) {
	print FINAL substr ($final_positive_contigs->{$name}, $position, 60), "\n";
    }
}
close FINAL;

print "----------------------------------------------\n";
print $log_file_handle "----------------------------------------------\n";

print_current_time("Generate positive contigs file. Check $output/final_result_dir/final_positive_contigs.fasta");

#
# if mapping, generate SAM file using bowtie2
#

if ($mapping =~ m/yes/) {
    print "Performing final mapping and generating SAM file using bowtie2... ";
    print_current_time("Performing final mapping and generating SAM file using bowtie2.");

    open(READS_LIST, ">$output/final_result_dir/reads_list.txt") or
        die "ERROR: Could not create file $output/final_result_dir/reads_list.txt (731): $!\n";
    foreach my $read (keys %$all_recruited_reads) {
        print READS_LIST $read, "\n";
    }
    close READS_LIST;

    my $bowtie_option_q_or_f;
    if($db_info->{'db_type'} =~ m/sff/i || $db_info->{'db_type'} =~ m/fasta/i) {
    	$bowtie_option_q_or_f = "-f";
    }

	my $cmd1 = "grep -F --no-group-separator -A 1 -f $output/final_result_dir/reads_list.txt $db > $output/final_result_dir/reads_list.$ext";
	`$cmd1`;
    !system "bowtie2-build $output/final_result_dir/final_positive_contigs.fasta $output/final_result_dir/final_positive_contigs.fasta 2>>$log >>$output/bowtie2-build.log"
	or warn "WARNING: bowtie2-build failed: $!\nCommand: bowtie2-build $output/final_result_dir/final_positive_contigs.fasta $output/final_result_dir/final_positive_contigs.fasta 2>>$log >>$output/bowtie2-build.log\n";
    !system "bowtie2 $bowtie_option_q_or_f -S $output/final_result_dir/final_positive_contigs.sam --fullref -a $output/final_result_dir/final_positive_contigs.fasta $output/final_result_dir/reads_list.$ext 2>>$output/bowtie2.log" or
	warn "WARNING: bowtie2 failed: $!\nCommand: bowtie2 $bowtie_option_q_or_f -S $output/final_result_dir/final_positive_contigs.sam --fullref -a $output/final_result_dir/final_positive_contigs.fasta $output/final_result_dir/reads_list.$ext 2>>$output/bowtie2.log\n";

    #
    # remove intermediate files
    #
    if ($clean =~ m/yes/i) {
	unlink "$output/bowtie2.log";
	unlink "$output/bowtie2-build.log";
	unlink glob "$output/final_result_dir/reads_list*";
	unlink glob "$output/final_result_dir/final_positive_contigs.fasta.*";
    }
    print "Done.\n";
}

#
# finish genseed-hmm and remove some intermediate files
#
unlink "$output/contigs_next_round.fasta" if (-e "$output/contigs_next_round.fasta");
unlink "$output/contigs.fa" if (-e "$output/contigs.fa");
unlink "$output/contigs.fasta" if (-e "$output/contigs.fasta");
unlink "$output/input_sequences.fasta.cap.contigs.qual" if (-e "$output/input_sequences.fasta.cap.contigs.qual");
unlink "$output/recruited_reads.fasta.cap.contigs.qual" if (-e "$output/recruited_reads.fasta.cap.contigs.qual");
unlink "$output/454AllContigs.qual" if (-e "$output/454AllContigs.qual");
unlink "$output/final_result_dir/final_contigs.fasta" if (-e "$output/final_result_dir/final_contigs.fasta" && $clean =~ m/yes/i);
unlink "$output/check_progress.txt" if ($clean =~ m/yes/i);
unlink "$output/makeblastdb.log" if (-e "$output/makeblastdb.log" && $clean =~ m/yes/i);
unlink glob "$output/*.qual" if (-e "$output/*.qual");
system "mv $output/check_progress.txt $output/check_progress_dir/final_check.txt" if (-e "$output/check_progress.txt" && $clean =~ m/no/i);
unlink glob "$output/contigs_assembler.fa.qual*" if (-e "$output/contigs_assembler.fa.qual");

my $duration = time - $start_time;
$duration = $duration/60;
$duration = sprintf("%.2f", $duration);
print_current_time("Execution time: $duration minutes");
print "\nExecution time: $duration minutes\n";
exit;

####################################
# subroutines

sub output_dir_name {
    my $output_dir_name = shift;
    my $count = 2;
    my $flag = 0;

    if (-d $output_dir_name) {
	$flag = 1;
	while (-d "$output_dir_name\_$count") {
	    $count++;
	}
	$output_dir_name = "$output_dir_name\_$count";
    }
    print "\nOutput directory already exists, saving analysis to $output_dir_name instead.\n\n" if $flag;

    return ($output_dir_name);
}

sub seed_type_verify {
    my $seed_name = shift;
    my $seed_type = 'DNA';

    open (DATA, $seed_name) or die "ERROR: Could not open seed file $seed_name (803): $!\n";

  TYPE_LOOP:
    while (my $line = <DATA>) {
	if ($line =~ m/^\s*$/) {
	    next TYPE_LOOP;
	}
	if ($line =~ m/^>/) {
	    next TYPE_LOOP;
	}

	#
	# HMM seed?
	#
	if ($line =~ m/^HMMER.*/) {
	    $seed_type = 'Profile HMM';
	    last TYPE_LOOP;
	}

	#
	# AA or DNA?
	#
	if ($line =~ m/([lvipfsyqderkhwm]+)/i) {
	    $seed_type = 'Protein';
	}
    }
    close DATA;
    return ($seed_type)
}

sub make_output_dir {
    system "mkdir $output";
    system "mkdir $output/final_result_dir";

    if ($clean =~ m/no/i) {
	system "mkdir $output/fasta_dir";
	system "mkdir $output/similarity_search_dir";
	system "mkdir $output/check_progress_dir";
	system "mkdir $output/assembler_dir";
    }
}

sub format_multifasta_headers {
    my $db = shift;
    my $given_name = shift;

    my $name;
    my $extension = "";

    if ($db =~ m/(.*)\.([\w]+)$/) {
        $name = $1;
        $extension = $2;
    }
    else {
	$name = $db;
    }

    my $formatted_headers_db = $name . "_formatted." . $extension;
    $formatted_headers_db = $given_name if defined $given_name;

    if (!-e $formatted_headers_db) {
	print "Eliminating duplicate headers from $db and creating $formatted_headers_db ... ";
	print_current_time("Eliminating duplicate headers from $db and creating $formatted_headers_db ... ");

	open(DDB, "<$db") or die "ERROR: Could not open database file $db (865): $!\n";
	open(FORMAT, ">$formatted_headers_db") or die "ERROR: Could not create file $formatted_headers_db (866): $!\n";

	my %seen;
	while(my $line = <DDB>) {
	    if($line =~ /^(>\S+)/) {
		my $id = $1;

		if(!exists $seen{$id}) {
		    print FORMAT "$id\_duplicate1\n";
		}
		else {
		    print FORMAT "$id\_duplicate2\n";
		}
		$seen{$id}++;
	    }
	    else {
		print FORMAT $line;
	    }
	}

	print "Done.\n";
	close DDB;
	close FORMAT;
    }
    else {
	print "File $formatted_headers_db already exists; skipping duplicate header removal.\n";
	print_current_time("File $formatted_headers_db already exists; skipping duplicate header removal.");
    }
    return ($formatted_headers_db);
}

sub run_makeblastdb {
    my $db = shift;

    if($db_info->{'db_type'} =~ m/fasta/i) {
	#
	# if BLAST-formatted DB does not exist, create it
	#
        unless (-e "$db\.nal" || -e "$db\.nhr") {
	    print "Generating BLAST index from $db ... ";
	    print_current_time("Generating BLAST index from $db ... ");
	    !system "makeblastdb -in $db -dbtype nucl -parse_seqids -logfile $output/makeblastdb.log"
		or die "ERROR: Formatting of $db by makeblastdb failed (897): $!\n";
	    print "Done.\n";
        }
	#
	# index DB file
	#
	($db_info->{'seq_cidx'}, $db_info->{'seq_cidxin'}, $id_table) = FA_index($db);

	#
	# if a FASTA quality file was given, index for retrieval
	#
        if($fasta_qual ne ' ') {
	    ($db_info->{'qual_cidx'}, $db_info->{'qual_cidxin'}, $qid_table) = FA_index($fasta_qual);
        }
	$db_info->{'blast_db'} = $db;
    }
    elsif($db_info->{'db_type'} eq 'sff') {
        #
        # there could be multiple SFF files in $db, therefore split $db and iterate
        #
        my @files = split(/\s+/, $db);
        unless(-e "SFFextracted.fasta.nal" || -e "SFFextracted.fasta.nhr") {
	    #
	    # if those files exist, database is already formatted, for example in a previous run;
	    # otherwise, create FASTA file with reads from all SFFs given and format for BLAST
	    #
	    foreach my $file (@files) {
		!system "sffinfo -s $file >> SFFextracted.fasta"
		    or die "ERROR: Could not convert $file to FASTA and append to SFFextracted.fasta (929): $!\n";
	    }
	    print "\nFormatting SFFextracted.fasta and generating BLAST index ... ";
	    print_current_time("Formatting SFFextracted.fasta and generating BLAST index ... ");
	    !system "makeblastdb -in SFFextracted.fasta -dbtype nucl -parse_seqids -logfile $output/makeblastdb.log"
		or die "ERROR: Formatting of SFFextracted.fasta by makeblastdb failed (933): $!\n";
	    print "Done.\n";
        }
        $db_info->{'blast_db'} = "SFFextracted.fasta";
    }
}

sub run_transeq {
    my $db_name = shift;
    my $transeq_db_name = shift;
    my $threads = shift;

    print "Performing conceptual translation of $db_name and creating $transeq_db_name... \n";
    print_current_time("Performing conceptual translation of $db_name and creating $transeq_db_name... ");

    if ($threads == 1) {
	system "transeq -sformat pearson $db_name $transeq_db_name -frame 6 2>>$log";
	print "Done.\n";
	return;
    }

    #
    # creating temporary directory
    #
    my $tmp_dir = $output . "/tmp_dir";
    mkdir "$tmp_dir"
        or die "ERROR: Could not create temporary directory $tmp_dir: $!\n";

    my $temporary_files_prefix = "$tmp_dir/run_transeq_temp_file_" . int(rand()*100) . "_";

    #
    # count records in multifasta database
    #
    my $fasta_quantity = `grep -c ">" $db_name`;
    my $limit = int($fasta_quantity/$threads)+ 1;

    local $/ = ">";

    open (IN, "<$db_name")
	or die "ERROR: Could not open file $db_name for run_transeq: $!\n";
    <IN>;

    my $i = 1;
    open (FILE,">$temporary_files_prefix" . "1.fasta")
	or die "ERROR: Could not create file $temporary_files_prefix" . "1.fasta: $!\n";
    my $num_seq = 0;

  PRINTLOOP:
    while (my $entry = <IN>) {
	chomp($entry);
	print FILE ">$entry";
	$num_seq++;
	if ($num_seq > $limit) {
	    close(FILE);
	    system "(transeq -sformat pearson $temporary_files_prefix$i.fasta $temporary_files_prefix$i.fasta.transeq -frame 6 2>>$log ; touch $tmp_dir/end_$i)&";
	    $i++;
	    open(FILE, ">$temporary_files_prefix$i.fasta")
		or die "ERROR: problem creating file $temporary_files_prefix$i.fasta: $!\n";
	    $num_seq = 0;
	}
    }

    if ($num_seq < $limit) {
	close(FILE);
	system "(transeq -sformat pearson $temporary_files_prefix$i.fasta $temporary_files_prefix$i.fasta.transeq -frame 6 2>>$log ; touch $tmp_dir/end_$i)&";
    }

    #
    # wait for background transeq runs and then concatenate all files in one large translated database
    #
    my $cnt = 0;
    print "Waiting for transeq jobs to finish... ";

    while($cnt < $i) {
	sleep(10);
	$cnt = 0;
	for(my $j=1; $j <= $i; $j++) {
	    if(-e "$tmp_dir/end_$j") {
		$cnt++;
	    }
	}
    }
    print "Done.\n";

    print "Concatenating temporary transeq output files... ";
    system "cat $temporary_files_prefix*.transeq > $transeq_db_name";
    system "rm -rf $tmp_dir";
    print "Done.\n";
}

sub print_current_time {
    my $message_for_print = shift;

    my @time = localtime(time);
    printf $log_file_handle "[%02d:%02d:%02d] $message_for_print\n", $time[2], $time[1], $time[0];
}

sub progressive_assembly_routine {
    my $last_round_contig_file = $seed_information->{'original_name'};

    my $round = 0;

    my $next_round_raw_contig_file;
    my $next_round_contig_file;
    my $next_round_contig_quality_file = " ";
    my $next_round_recruited_reads;
    my $all_recruited_reads = {};
    my $current_contigs = {};
    my $stop_progressive_assembly = 0;

  PROGRESSIVE_ASSEMBLY_LOOP:
    while (1) {
	$round++;

	print "\n###     Round $round     ###\n";
	print $log_file_handle "###     Round $round     ###\n";

	if ($max_number_rounds) {
	    if ($round > $max_number_rounds) {
		print $log_file_handle "\nGenSeed-HMM reached the maximum number of rounds. \n";
		print "GenSeed-HMM reached the maximum number of rounds.\n";
		unlink "$output/recruited_reads.fasta.cap.contigs.qual" if (-e "$output/recruited_reads.fasta.cap.contigs.qual");
		last PROGRESSIVE_ASSEMBLY_LOOP;
	    }
	}
	$next_round_recruited_reads = "$output/recruited_reads.$ext";

	#
	# first recruit new reads
	#
	my $new_reads_in_this_round;

	$stop_progressive_assembly = recruit_new_reads(\$new_reads_in_this_round,
						       $last_round_contig_file, 
						       $next_round_recruited_reads, 
						       $all_recruited_reads, 
						       $round, 
						       $next_round_contig_quality_file);
	#
	# fix recruited_reads sequence IDs (remove lcl| -- or whatever 3 letters appear after > and before |)
	#
	system("sed -i 's/>...\|/>/' $next_round_recruited_reads") if (-e $next_round_recruited_reads);

	#
	# check if there are new reads in this round and store the intermediate files if user wants
	#
	if ($clean =~ m/no/i) {
	    system "mv $output/tabular_temp.results $output/similarity_search_dir/result_similarity_search_round$round.txt" if(-e "$output/tabular_temp.results");
	    system "mv $output/seed_file.fasta $output/fasta_dir/seed_file_round$round.fasta" if (-e "$output/seed_file.fasta");
	    unlink "$output/tabular.results" if (-e "$output/tabular.results");
	}

	if ($stop_progressive_assembly) {
	    open(FINAL_FILE, ">>$output/final_result_dir/final_contigs.fasta") or die "ERROR: Could not write to file $output/final_result_dir/final_contigs.fasta (1109): $!\n";
	    foreach my $key (keys %$current_contigs) {
		print "Contig $key did not recruit reads - finishing its progressive assembly\n";
		print_current_time("Contig $key did not recruit reads - finishing its progressive assembly");
	    }
	    close FINAL_FILE;

	    last PROGRESSIVE_ASSEMBLY_LOOP;
	}

	#
	# dereplicate $next_round_recruited_reads
	#
	if ($assembler_name =~ m/cap3/i) {
		print "There is $new_reads_in_this_round new reads. Dereplicating the reads\n";
		print_current_time("There is $new_reads_in_this_round new reads. Dereplicating the reads");
		my $dereplicated_reads = "$output/dereplicated_recruited_reads.$ext";
		my $cmdDerep = "vsearch --derep_fulllength $next_round_recruited_reads --notrunclabels --output $dereplicated_reads --quiet";
		`$cmdDerep`;
		system "mv $dereplicated_reads $next_round_recruited_reads";
		my $new_count = `grep -c ">" $next_round_recruited_reads`;
		chomp($new_count);
		$new_reads_in_this_round = $new_count + 0;
	}
	#
	# now run the progressive assembly round
	#
	print "Running progressive assembly with: $new_reads_in_this_round dereplicated new reads\n";
	print_current_time("Running progressive assembly with: $new_reads_in_this_round dereplicated new reads");

	if ($assembler_name =~ m/newbler/i) {
	    $stop_progressive_assembly = run_newbler($last_round_contig_file, 
						     $next_round_recruited_reads, 
						     \$next_round_raw_contig_file, 
						     \$next_round_contig_quality_file, 
						     $round);
	}
	if ($assembler_name =~ m/cap3/i) {
	    $stop_progressive_assembly = run_cap3($last_round_contig_file, 
						  $next_round_recruited_reads, 
						  \$next_round_raw_contig_file, 
						  \$next_round_contig_quality_file);
	}
	if ($assembler_name =~ m/abyss/i) {
	    $stop_progressive_assembly = run_abyss($last_round_contig_file, 
						   $next_round_recruited_reads, 
						   \$next_round_raw_contig_file, 
						   $round);
	}
	if ($assembler_name =~ m/velvet/i) {
	    $stop_progressive_assembly = run_velvet($last_round_contig_file, 
						    $next_round_recruited_reads, 
						    \$next_round_raw_contig_file, 
						    $round);
	}
	if ($assembler_name =~ m/soapdenovo/i) {
	    $stop_progressive_assembly = run_soapdenovo($last_round_contig_file, 
							$next_round_recruited_reads, 
							\$next_round_raw_contig_file, 
							$round);
	}

	#
	# checks if assembly changed after this round and store intermediate files if user wants and remove intermediate files
	#
	if ($seed_information->{'type'} !~ m/DNA/ && $clean =~ m/no/) {
	    system "mv $output/recruited_reads.$ext $output/assembler_dir/input_for_assembler_round$round.$ext";
	    system "mv $output/recruited_reads.$ext\.qual $output/assembler_dir/input_for_assembler_round$round.$ext\.qual" if (-e "$output/recruited_reads.$ext\.qual");
	}
	else {
	    unlink "$output/recruited_reads.$ext" if (-e "$output/recruited_reads.$ext");
	    unlink "$output/recruited_reads.$ext\.qual" if (-e "$output/recruited_reads.$ext\.qual");
	}

	if ($stop_progressive_assembly) {
	    if ($round == 1) {
		print "Assembler could not reconstruct any contig(s). Aborting genseed-hmm.pl\n";
		exit;
	    }
	    else {
		print "Assembler could not reconstruct any contig(s).\n";
		last PROGRESSIVE_ASSEMBLY_LOOP;
	    }
	}

	if (-e $next_round_raw_contig_file && $next_round_raw_contig_file) {
	    system "mv $next_round_raw_contig_file $output/contigs_assembler.fa";
	    $next_round_raw_contig_file = "$output/contigs_assembler.fa";
	}
	if (-e $next_round_contig_quality_file && $next_round_contig_quality_file) {
	    system "mv $next_round_contig_quality_file $output/contigs_assembler.fa.qual";
	    $next_round_contig_quality_file = "$output/contigs_assembler.fa.qual";
	    #
	    # fix quality identifiers to match those from sequence file
	    #
	    system "sed -i 's/>/>r$round./' $next_round_contig_quality_file";
	}

	my $positive_contigs = {};
	$stop_progressive_assembly = check_progress($last_round_contig_file, 
						    $next_round_raw_contig_file, 
						    $seed_information->{'type'}, 
						    $output, 
						    $positive_contigs,
						    $current_contigs,
						    $round);

	if ($clean =~ m/no/) {
	    my $last_round = $round-1;
	    system "mv $next_round_raw_contig_file $output/assembler_dir/reconstructed_contig_round$round.fasta" if (-e $next_round_raw_contig_file);
	    system "cp $next_round_contig_quality_file $output/assembler_dir/reconstructed_contig_round$round.fasta.qual" if (-e $next_round_contig_quality_file);
	    system "mv $output/input_sequences.$ext $output/assembler_dir/input_for_assembler_round$round.$ext" if (-e "$output/input_sequences.$ext");
	    system "mv $output/input_sequences.$ext.qual $output/assembler_dir/input_for_assembler_round$round.$ext.qual" if (-e "$output/input_sequences.$ext.qual");
	    system "mv $output/assembler.log $output/assembler_dir/assembler_log_round$round.txt" if (-e "$output/assembler.log");
	    system "mv $output/check_progress.txt $output/check_progress_dir/check_contig_progress_round$round.txt" if (-e "$output/check_progress.txt");
	    system "mv $next_round_contig_file $output/fasta_dir/positive_contigs_round$last_round.fasta" if ($next_round_contig_file);
	}
	else {
	    unlink "$next_round_raw_contig_file" if (-e $next_round_raw_contig_file);
	    unlink "$next_round_contig_quality_file" if (-e $next_round_contig_quality_file);
	}

	if ($stop_progressive_assembly) {
	    last PROGRESSIVE_ASSEMBLY_LOOP;
	}

	#
	# print the positive contigs to a file
	#
	$next_round_contig_file = "$output/contigs_next_round.fasta";

	open(POSITIVE, ">$next_round_contig_file") or die "ERROR: Could not create file $next_round_contig_file (1303): $!\n";

	foreach my $key (sort { length $positive_contigs->{$a} <=> length $positive_contigs->{$b} } keys %$positive_contigs) {
	    my $length = length $positive_contigs->{$key};
	    print "Length of positive contig $key: $length bp\n";
	    print_current_time("Reconstructed positive contig $key with: $length bp");
	    print POSITIVE ">$key\n";
	    print POSITIVE $positive_contigs->{$key}, "\n";
	}
	close POSITIVE;

	#
	# if a contig for whatever reason does not get used in an assembly (i.e. ends up as a singleton), it gets silently dropped from the analysis;
	# the following is intended to detect such contigs (by checking whether old contigs are similar to new contigs -- if they are not, then save them)
	#
	if(keys %$current_contigs) {
	    my(@fail) = compare_oldnew_contigs($current_contigs, $positive_contigs) if $round > 1;

	    if(@fail) {
		foreach my $ctg (@fail) {
		    open(FINAL_FILE, ">>$output/final_result_dir/final_contigs.fasta") or die "ERROR: Could not write to file $output/final_result_dir/final_contigs.fasta (1224): $!\n";
		    print "Contig $ctg did not return from assembly - finishing its progressive assembly\n";
		    print_current_time("Contig $ctg did not return from assembly - finishing its progressive assembly");
		    print FINAL_FILE ">$ctg\n$current_contigs->{$ctg}\n";
		    close FINAL_FILE;
		}
	    }
	}

	#
	# change seed type to DNA every round (in case the seed is an HMM-profile or protein)
	# store current contigs in a hash
	#
	$current_contigs = $positive_contigs;
	$last_round_contig_file = $next_round_contig_file;
	$seed_information->{'type'} = 'DNA';
    }
    return($all_recruited_reads, $current_contigs);
}

sub compare_oldnew_contigs {
    my $old_contigs = shift;
    my $new_contigs = shift;

    open(TMPN, ">$output/cmp_ctg_N") or die "ERROR: Could not create temporary file $output/cmp_ctg_N (1243): $!\n";
    open(TMPO, ">$output/cmp_ctg_O") or die "ERROR: Could not create temporary file $output/cmp_ctg_O (1244): $!\n";
    foreach my $key (keys %$new_contigs) {
	print TMPN ">$key\n$new_contigs->{$key}\n";
    }
    foreach my $key (keys %$old_contigs) {
	print TMPO ">$key\n$old_contigs->{$key}\n";
    }
    close TMPN;
    close TMPO;

    !system("blastn -query $output/cmp_ctg_N -subject $output/cmp_ctg_O -outfmt 6 -out $output/cmp_out -perc_identity 95") or
	die "ERROR: Could not run BLAST for contig comparison (1255): $!\n";

    open(CMP, "<$output/cmp_out") or die "ERROR: Could not open file $output/cmp_out (1257): $!\n";
    my %subj;

    while(my $line = <CMP>) {
	next if $line =~ /^\s*$/;
	my(@data) = split(/\t/, $line);
	$subj{$data[1]}++;
    }

    my @fail;
    foreach my $res (keys %$old_contigs) {
	unless(exists $subj{$res}) {
	    push @fail, $res;
	}
    }

    unlink "$output/cmp_ctg_N";
    unlink "$output/cmp_ctg_O";
    unlink "$output/cmp_out";

    return @fail;
}

sub recruit_new_reads {
    my $new_reads_in_this_round = shift;
    my $sequence_for_similarity_search = shift;
    my $recruited_reads_filename = shift;
    my $all_recruited_reads = shift;
    my $round = shift;
    my $contig_quality_file = shift;

    my $sequence_type = $seed_information->{'type'};

    #
    # if round > 1, generate the seed file by extracting the contig(s) ends 
    #
    if ($round > 1) {
	print_current_time("Generating seed file for similarity search");
	my $next_seed_file = "$output/seed_file.fasta";
	generate_sequence_for_similarity_search($sequence_for_similarity_search, $next_seed_file, $contig_quality_file);
	$sequence_for_similarity_search = $next_seed_file;
    }

    #
    # run similarity search with hmmsearch (seed = Profile HMM), blastn (seed = DNA) or tblastn (seed = Protein)
    # check if hits exist and, if so, retrieve their sequences
    #
    if ($sequence_type =~ m/Profile HMM/) {
	print_current_time("Performing similarity search");
	!system "hmmsearch $hmmsearch_parameters --tblout $output/tabular_temp.results -o $output/full_temp.results $sequence_for_similarity_search $db_info->{'transeq'} 2>> $log" or
	    die "ERROR: Could not run hmmsearch (1237): $!\nCommand: hmmsearch $hmmsearch_parameters --tblout $output/tabular_temp.results -o $output/full_temp.results $sequence_for_similarity_search $db_info->{'transeq'} 2>> $log\n";
	unlink "$output/full_temp.results";
    }
    elsif ($sequence_type =~ m/DNA/) {
	print_current_time("Performing similarity search");
	!system "blastn -query $sequence_for_similarity_search -db $db_info->{'blast_db'} $blastn_parameters -outfmt \"6 sseqid qlen length\" -out $output/tabular_temp.results 2>>$log" or
	    die "ERROR: Could not run blastn (1243): $!\nCommand: blastn -query $sequence_for_similarity_search -db $db_info->{'blast_db'} $blastn_parameters -outfmt \"6 sseqid qlen length\" -out $output/tabular_temp.results 2>>$log\n";
    }
    else {
	print_current_time("Performing similarity search");
	!system "tblastn -query $sequence_for_similarity_search -db $db_info->{'blast_db'} $tblastn_parameters -outfmt \"6 sseqid qlen length\" -out $output/tabular_temp.results 2>>$log" or
	    die "ERROR: Could not run tblastn (1248): $!\nCommand: tblastn -query $sequence_for_similarity_search -db $db_info->{'blast_db'} $tblastn_parameters -outfmt \"6 sseqid qlen length\" -out $output/tabular_temp.results 2>>$log\n";
    }

    print_current_time("Filtering hits from similarity search");

    open(SIMILARITY_FILE, "$output/tabular_temp.results") or die "ERROR: Could not open file $output/tabular_temp.results (1253): $!\n";

    my @recruited_reads = ();

    if ($sequence_type =~ /Profile HMM/i) {
	while(my $line = <SIMILARITY_FILE>) {
	    if ($line =~ m/(\S+)(_\d)\s/) {
		push (@recruited_reads, $1);
	    }
	    
	}
    }
    else {
	while (my $line = <SIMILARITY_FILE>) {
	    my @filter_blast = split(/\t/, $line);
	    if ((($filter_blast[2]/$filter_blast[1])*100) >= $min_hsp_length) {
		push (@recruited_reads, $filter_blast[0]);
	    }
	}
    }
    close SIMILARITY_FILE;

    if (!@recruited_reads) {
	if ($round == 1) {
	    print "No match found with selected database! Select another seed, check the parameters or change them in similarity search (blastn, tblastn or hmmsearch).\n";
	    exit;
	}
	else {
	    print "No new read incorporated by GenSeed-HMM.\n";
	    return(1);
	}
    }

    open (FILE, ">$output/tabular.results") or die "ERROR: Could not create file $output/tabular.results (1286): $!\n";

    my $stop_progressive_assembly = 1;

    for my $hit_name (@recruited_reads) {
	if (!exists $all_recruited_reads->{$hit_name}) {
            $all_recruited_reads->{$hit_name}++;
            $$new_reads_in_this_round++;
            $stop_progressive_assembly = 0;
	    #
	    # redundancy can be avoided without comprimising the assembly process
	    #
	    if ($assembler_name =~ m/cap3/) {
		print FILE $hit_name, "\n";
	    }
        }
    }

    #
    # redundancy and coverage are important for Newbler, ABySS, SOAPdenovo and Velvet
    #
    if ($assembler_name !~ m/cap3/i) {
	foreach my $hit_name (@recruited_reads) {
	    print FILE $hit_name, "\n";
	}
    }
    close FILE;

    if ($stop_progressive_assembly) {
	print_current_time("No new read incorporated by GenSeed-HMM. ");
	print "No new read incorporated by GenSeed-HMM. \n";
	return (1);
    }
    else {

	#
	# retrieve new sequences 
	#
	print_current_time("Retrieving new reads for progressive assembly");
	my $cmd1 = "grep -F --no-group-separator -A 1 -f $output/tabular.results $db > $recruited_reads_filename";
	`$cmd1`;
    }

    #
    # remove intermediate files
    #
    system "rm $output/tabular.results";
}

sub generate_sequence_for_similarity_search {
    my $sequence = shift;
    my $next_seed_file = shift;
    my $contig_quality_file = shift;

    my %start_end_seed;

    if (($contig_quality_file and !-e $contig_quality_file)) {
	$contig_quality_file = undef;
    }

    if ($contig_quality_file and @quality_inspector) {
	#
	# extract contig ends, relying in contig qualities (if -quality_inspector is being used), and create the extension seed sequence for similarity search
	#
	my $length_quality_inspector = $quality_inspector[0];
	my $minimum_quality_value = $quality_inspector[1];
	my $maximum_percentage_bad_quality = $quality_inspector[2];

	my %quality_contigs;
	my $contig_name;

	open(QUALITY_FILE, $contig_quality_file) or die "ERROR: Could not open quality file $contig_quality_file (1357): $!\n";

	while (my $line = <QUALITY_FILE>) {
	    if ($line =~ m/^(>.+)/) {
		chomp $1;
		$contig_name = $1;
	    }
	    else {
		my @temporary_array = split (/\s/, $line);
		push(@{$quality_contigs{$contig_name}}, @temporary_array);
	    }
	}

	#
	# Start quality inspector
	#
	my $end_position;

	$length_quality_inspector--;

	for my $name (sort keys %quality_contigs) {
	    my $five_line_position = -1;
            my $three_line_position = -1;

	  QUALITY_LOOP:
	    for (my $start_position = 0; $start_position < scalar (@{$quality_contigs{$name}} - $length_quality_inspector); $start_position++) {
		$end_position = $start_position + $length_quality_inspector;

		my $bad_quality_value = 0;

		my @sliding_window = @{$quality_contigs{$name}}[$start_position..$end_position];

		#
		# check for bases with quality below the threshold
		#
		for my $quality_value (@sliding_window) {
		    if ($quality_value < $minimum_quality_value) {
			$bad_quality_value++;
		    }
		}

		if ( ($bad_quality_value/$length_quality_inspector) > $maximum_percentage_bad_quality/100) {
		    next QUALITY_LOOP;
		}

		#
		# create the start and end positions based on quality values
		#
		$three_line_position = $start_position + $length_quality_inspector;

		if ($five_line_position == -1) {
		    $five_line_position = $start_position;
		}
	    }

	    #
	    # store this information in a hash
	    #
	    $start_end_seed{$name}{'start'} = $five_line_position;
	    $start_end_seed{$name}{'end'} = $three_line_position+1;
	}
    }

    #
    # extract contig ends and create the extension seed sequence for similarity search
    #
    open(CONTIGS, $sequence) or die "ERROR: Could not open contig file $sequence (1423): $!\n";
    open(SEED, ">$next_seed_file") or die "ERROR: Could not create extension seed file $next_seed_file (1424): $!\n";
    
    my $header_name;
    
    while(my $line = <CONTIGS>) {
	if ($line =~ m/>/) {
	    chomp $line;
	    $header_name = $line;
	}
	else {
	    chomp $line;

	    #
	    # if extension seed size is bigger than contig, use entire contig as a seed in this round
	    #
	    if ($next_seed_size > length $line) {
		print SEED "$header_name\_full_seed_sequence\n";
		print SEED $line, "\n";
	    }
	    else {

		#
		# if extension seed size is smaller than contig, extract the contig ends and use these ends as seeds in this round
		#
		if (!%start_end_seed) {
		    if ($exp_direction =~ m/left|both/i) {
			print SEED "$header_name\_five_line_end\n";
			print SEED substr ($line, 0, $next_seed_size), "\n";
		    }
		    if ($exp_direction =~ m/right|both/i) {
			print SEED "$header_name\_three_line_end\n";
			print SEED substr ($line, ( (length $line) - $next_seed_size ), length $line), "\n";
		    }
		}
		else {
		    foreach my $name (keys %start_end_seed) {
			#
			# create the seed file based on the quality values
			#
			if ($name =~ m/$header_name/) {
			    if ($exp_direction =~ m/left|both/i) {
				print SEED "$header_name\_five_line_end\n";
				print SEED substr ($line, $start_end_seed{$name}{'start'}, $next_seed_size), "\n";
			    }
			    if ($exp_direction =~ m/right|both/i) {
				print SEED "$header_name\_three_line_end\n";
				print SEED substr ($line, ( ($start_end_seed{$name}{'end'}) - $next_seed_size ), $start_end_seed{$name}{'end'}), "\n";
			    }
			}
		    }
		}
	    }
	}
    }
    close CONTIGS;
    close SEED;
}

sub run_newbler {
    my $last_round_contig_file = shift;
    my $recruited_reads = shift;
    my $next_round_raw_contig_file = shift;
    my $next_round_contig_quality_file = shift;
    my $round = shift;

    my $input_file_for_assembler = " ";

    $input_file_for_assembler = get_file($recruited_reads, $last_round_contig_file, "newbler", $next_round_contig_quality_file, $round);
    if($db_info->{'db_type'} ne 'sff') { $recruited_reads = ' '; }

    unlink "$input_file_for_assembler\.qual" unless $use_qual =~ /yes/i;

    #
    # Newbler demands that reads and quality records be in the same order; therefore, here we order reads to be in the same order as the quality values
    #
    if(-e "$input_file_for_assembler\.qual") {
	open(TMPSORT, ">$output/tmp_sort") or die "ERROR: Could not create temporary file $output/tmp_sort (1669): $!\n";
	my(@qual_ids) = `grep \">\" $input_file_for_assembler\.qual`;
	my %empty_list;

	foreach my $id (@qual_ids) {
	    chomp $id;
	    my(@id) = split(/\s+/, $id);
	    $id[0] = substr($id[0], 1);
	    my $seq = return_sequence_from_fasta($input_file_for_assembler, $id[0]);
	    unless($seq eq "-1") {
		print TMPSORT ">$id[0]\n$seq\n";
	    }
	    else {
		$empty_list{$id[0]}++;
	    }
	}
	close TMPSORT;

	#
	# remove qual scores for seqs that did not make to here (saved in %empty_list)
	#
	open(Q, "<$input_file_for_assembler\.qual") or die "ERROR: Could not open file $input_file_for_assembler\.qual (1688): $!\n";
	open(TMPQ, ">$output/tmp_q") or die "ERROR: Could not create temporary file $output/tmp_q (1689): $!\n";
	my $qflag = 0;
	while(<Q>) {
	    if($_ =~ /^>/) {
		chomp;
		my(@id) = split(/\s+/, $_);
		$id[0] = substr($id[0], 1);
		if(exists $empty_list{$id[0]}) {
		    $qflag = 0;
		}
		else {
		    $qflag = 1;
		}
	    }
	    print TMPQ "$_\n" if $qflag;
	}
	close TMPQ;
	system "mv $output/tmp_sort $input_file_for_assembler";
	system "mv $output/tmp_q $input_file_for_assembler\.qual";
    }

    #
    # run newbler
    #
    !system "runAssembly -m -urt -o $output/newbler $assembler_parameters $recruited_reads $input_file_for_assembler >>$output/assembler.log 2>>$log"
	or die "ERROR: Newbler assembly died: $!\nCheck assembler log file.\n";

    remove_intermediary_assembler_files();

    my $stop_progressive_assembly = check_assembler_reconstruction("$output/454AllContigs.fna");

    if ($stop_progressive_assembly) {
	return(1);
    }
    else {
	$$next_round_raw_contig_file = "$output/454AllContigs.fna";
	$$next_round_contig_quality_file = "$output/454AllContigs.qual";
	return(0);
    }
}

sub run_cap3 {
    my $last_round_contig_file = shift;
    my $recruited_reads= shift;
    my $next_round_raw_contig_file = shift;
    my $next_round_contig_quality_file = shift;

    my $input_file_for_assembler;

    $input_file_for_assembler = get_file($recruited_reads, $last_round_contig_file, "cap3");

    unlink "$input_file_for_assembler\.qual" unless $use_qual =~ /yes/i;

    #
    # run cap3
    #
    system "cap3 $input_file_for_assembler $assembler_parameters >$output/assembler.log 2>>$log";
    remove_intermediary_assembler_files($input_file_for_assembler);

    my $stop_progressive_assembly = check_assembler_reconstruction("$input_file_for_assembler.cap.contigs");

    if ($stop_progressive_assembly) {
	return(1);
    }
    else {
	$$next_round_raw_contig_file = "$input_file_for_assembler.cap.contigs";
	$$next_round_contig_quality_file = "$input_file_for_assembler.cap.contigs.qual";
	return(0);
    }
}

sub run_abyss {
    my $last_round_contig_file = shift;
    my $recruited_reads= shift;
    my $next_round_raw_contig_file = shift;
    my $round = shift;

    my $input_file_for_assembler;

    $input_file_for_assembler = get_file($recruited_reads, $last_round_contig_file, "abyss");

    #
    # run abyss
    #
    system "abyss -k $kmer $assembler_parameters -o $output/contigs.fasta $input_file_for_assembler >>$output/assembler.log 2>>$output/assembler.log";

    #
    # check the assembler reconstruction
    #
    my $stop_progressive_assembly = check_assembler_reconstruction("$output/contigs.fasta");

    if ($stop_progressive_assembly) {
	if ($round == 1) {

	    #
	    # force progressive assembly process decreasing the kmer and if necessary trimming the reads end
	    #
	    $stop_progressive_assembly = force_progressive_assembly($input_file_for_assembler);
	    
	    if ($stop_progressive_assembly) {
		print "ABySS cannot reconstruct any contig(s) in the first round. Aborting.\n";
		exit;
	    }
	    else {
		$$next_round_raw_contig_file = "$output/contigs.fasta";
		return(0);
	    }
	}
	else {
	    return(1);
	}
    }
    else {
	$$next_round_raw_contig_file = "$output/contigs.fasta";
	return(0);
    }
}

sub run_velvet {
    my $last_round_contig_file = shift;
    my $recruited_reads= shift;
    my $next_round_raw_contig_file = shift;
    my $round = shift;

    my $input_file_for_assembler;

    $input_file_for_assembler = get_file($recruited_reads, $last_round_contig_file, "velvet");

    #
    # run velvetg and velveth 
    #
    system "velveth $output/ $kmer $input_file_for_assembler >> $output/assembler.log 2>> $log";
    system "velvetg $output/ $assembler_parameters >> $output/assembler.log 2>> $log";

    #
    # remove intermediate velvet files
    #
    remove_intermediary_assembler_files();

    my $stop_progressive_assembly = check_assembler_reconstruction("$output/contigs.fa");

    if ($stop_progressive_assembly) {
        if ($round == 1) {

	    #
            # force progressive assembly process decreasing the kmer and if necessary trimming the reads end
	    #
	    $stop_progressive_assembly = force_progressive_assembly($input_file_for_assembler);

	    if ($stop_progressive_assembly) {
		print "Velvet cannot reconstruct any contig(s) in the first round. Aborting.\n";
		exit;
	    }
	    else {
		$$next_round_raw_contig_file = "$output/contigs.fa";
		return(0);
	    }
	}
        else {
	    return(1);
        }
    }
    else {
	$$next_round_raw_contig_file = "$output/contigs.fa";
	return(0);
    }
}

sub run_soapdenovo {
    my $last_round_contig_file = shift;
    my $recruited_reads= shift;
    my $next_round_raw_contig_file = shift;
    my $round = shift;

    my $input_file_for_assembler;
    my $soap_assembler;

    $input_file_for_assembler = get_file($recruited_reads, $last_round_contig_file, "soapdenovo");

    #
    # Create soapdenovo configurator file for this round
    #
    open(SOAP_CONF, ">$output/soap_configuration_file.temp") or die "ERROR: Could not create SOAP configuration file $output/soap_configuration_file.temp : $!\n";
    print SOAP_CONF "max_rd_len=1000000\n";
    print SOAP_CONF "[LIB]\n";
    print SOAP_CONF "asm_flags=1\n";
    print SOAP_CONF "f=$input_file_for_assembler";
    close SOAP_CONF;

    #
    # Run soapdenovo assembler. Support large kmer up to 127 to utilize long reads. Three versions are provided. (SOAPdenovo-127kmer, SOAPdenovo-63kmer, SOAPdenovo-31kmer)
    #
    if ($kmer <= 127 && $kmer > 63) {
	$soap_assembler = "SOAPdenovo-127mer";
    }
    elsif ($kmer <= 63 && $kmer > 31) {
	$soap_assembler = "SOAPdenovo-63mer";
    }
    elsif ($kmer <= 31) {
	$soap_assembler = "SOAPdenovo-31mer";
    }

    system "$soap_assembler pregraph -K $kmer -s $output/soap_configuration_file.temp $assembler_parameters -o $output/contigs >> $output/assembler.log 2>> $log";
    system "$soap_assembler contig -g $output/contigs >> $output/assembler.log 2>>$log";

    remove_intermediary_assembler_files();

    my $stop_progressive_assembly = check_assembler_reconstruction("$output/contigs.contig");

    if ($stop_progressive_assembly) {
        if ($round == 1) {
	    #
	    # force progressive assembly process decreasing the kmer and if necessary trimming the reads end
	    #
	    $stop_progressive_assembly = force_progressive_assembly($input_file_for_assembler);
	    
	    if ($stop_progressive_assembly) {
		print "SOAPdenovo cannot reconstruct any contig(s) in the first round. Aborting.\n";
		exit;
	    }
	    else {
		$$next_round_raw_contig_file = "$output/contigs.contig";
		return(0);
	    }
	}
        else {
	    return(1);
        }
    }
    else {
	$$next_round_raw_contig_file = "$output/contigs.contig";
	return(0);
    }
}

sub generate_newbler_input_file {
    my $last_round_contig_file = shift;
    my $input_sequences = shift;
    my $next_round_qual = shift;

    #
    # Newbler ignores input sequences that are longer than 1999 bp and does not retain sequences represented by a single read. 
    # In the context of progressive assembly process, this is a problem.
    # Based on the work of Laurent Keller et. al. 2011, we developed a similar approuch to solve this problem:
    # We split the sequences >1999 bp into subsequences of 1500 bp with a 1300 overlap using EMBOSS' splitter, 
    # and these sequences can be used by Newbler and GenSeed-HMM in the progressive assembly process.
    #
    my $ncf = "$output/newbler_contig.fasta";

    system "splitter $last_round_contig_file $ncf\.split -size 1500 -overlap 1300  2>>$log";
    
    for(my $i=0; $i < $redundancy; $i++) {
        system "cat $ncf\.split >> $input_sequences";
    }

    #
    # add Newbler contig quals to retrieved read quals, if user has given a FASTA quality file
    #
    if(defined $$next_round_qual) {
	if (-e $$next_round_qual && $fasta_qual ne ' ') {
	    my $osq = quality_splitter($$next_round_qual);
	    for(my $i=0; $i < $redundancy; $i++) {
		system "cat $osq >> $input_sequences\.qual"; 
	    }
	}
    }

    #
    # remove intermediate files
    #
    unlink glob "$ncf\*";
}

sub remove_intermediary_assembler_files {
    my $contig_cap3 = shift;

    if ($assembler_name =~ m/Newbler/i) {
	system "mv $output/newbler/454AllContigs.fna $output" if (-e "$output/newbler/454AllContigs.fna");
	system "mv $output/newbler/454AllContigs.qual $output" if (-e "$output/newbler/454AllContigs.qual");
	system "rm -rf $output/newbler";
    }
    elsif ($assembler_name =~ m/cap3/i) {
	system "rm $contig_cap3.cap.ace" if  (-e "$contig_cap3.cap.ace");
	system "rm $contig_cap3.cap.contigs.links" if (-e "$contig_cap3.cap.contigs.links");
	system "rm $contig_cap3.cap.info" if (-e "$contig_cap3.cap.info");
	system "rm $contig_cap3.cap.singlets" if (-e "$contig_cap3.cap.singlets");
    }
    elsif ($assembler_name =~ m/Velvet/i) {
	system "rm $output/Graph" if (-e "$output/Graph");
	system "rm $output/Graph2" if (-e "$output/Graph2");
	system "rm $output/LastGraph" if (-e "$output/LastGraph");
	system "rm $output/PreGraph" if (-e "$output/PreGraph");
	system "rm $output/Roadmaps" if (-e "$output/Roadmaps");
	system "rm $output/Log" if (-e "$output/Log");
	system "rm $output/stats.txt" if (-e "$output/stats.txt");
	system "rm $output/Sequences" if (-e "$output/Sequences");
    }
    elsif ($assembler_name =~ m/SOAPdenovo/i) {
	system "rm $output/soap_configuration_file.temp" if (-e "$output/soap_configuration_file.temp");
	system "rm $output/contigs.Arc" if (-e "$output/contigs.Arc");
	system "rm $output/contigs.ContigIndex" if (-e "$output/contigs.ContigIndex");
	system "rm $output/contigs.edge" if (-e "$output/contigs.edge");
	system "rm $output/contigs.kmerFreq" if (-e "$output/contigs.kmerFreq");
	system "rm $output/contigs.preArc" if (-e "$output/contigs.preArc");
	system "rm $output/contigs.preGraphBasic" if (-e "$output/contigs.preGraphBasic");
	system "rm $output/contigs.updated.edge" if (-e "$output/contigs.updated.edge");
	system "rm $output/contigs.vertex" if (-e "$output/contigs.vertex");
    }
}

sub check_assembler_reconstruction {
    my $reconstruction_file  = shift;

    #
    # Check if there are any contigs in the contig file 
    #
    if ((!-e "$reconstruction_file" || -z "$reconstruction_file")) {
	return(1);
    }
}

sub force_progressive_assembly {
    my $input_file_for_assembler = shift;

    #
    # first attempt: only decrease k-mer size and rerun the assembly
    #
    print_current_time("First attempt to initiate the progressive assembly: decreasing the k-mer size and rerun the assembly");
    print "First attempt to initiate the progressive assembly: decreasing k-mer size and rerun the assembly\n";

    my @decreasing_values = (2, 4, 6, 8);
    my $stop_progressive_assembly;

  FIRST_ATTEMPT:
    foreach my $value (@decreasing_values) {
	my $decreased_kmer = $kmer - $value;

	if ($assembler_name =~ m/ABySS/i) {
	    system "abyss -k $decreased_kmer $assembler_parameters -o $output/contigs.fasta $input_file_for_assembler >>$output/assembler.log 2>>$output/assembler.log";
	    $stop_progressive_assembly = check_assembler_reconstruction("$output/contigs.fasta");
	}
	elsif ($assembler_name =~ m/Velvet/i) {
	    system "velveth $output/ $decreased_kmer $input_file_for_assembler >> $output/assembler.log 2>> $log";
	    system "velvetg $output/ $assembler_parameters >> $output/assembler.log 2>> $log";

	    remove_intermediary_assembler_files();
	    $stop_progressive_assembly = check_assembler_reconstruction("$output/contigs.fa");
	}
	elsif ($assembler_name =~ m/SOAPdenovo/i) {
	    my $soap_assembler;
	    
	    if ($kmer <= 127 && $kmer > 63) {
		$soap_assembler = "SOAPdenovo-127mer";
	    }
	    elsif ($kmer <= 63 && $kmer > 31) {
		$soap_assembler = "SOAPdenovo-63mer"; 
	    }
	    elsif ($kmer <= 31) {
		$soap_assembler = "SOAPdenovo-31mer";
	    }
	    system "$soap_assembler pregraph -K $decreased_kmer -s $output/soap_configuration_file.temp $assembler_parameters -o $output/contigs >> $output/assembler.log 2>> $log";
	    system "$soap_assembler contig -g $output/contigs >> $output/assembler.log 2>>$log";
	    
	    remove_intermediary_assembler_files();
	    $stop_progressive_assembly = check_assembler_reconstruction("$output/contigs.contig");
	}
	if(!$stop_progressive_assembly) {
	    last FIRST_ATTEMPT;
	}
    }

    #
    # second attempt: decrease k-mer size and trim read ends (perhaps we have many bad quality bases in the end of reads)
    #
    if ($stop_progressive_assembly) {
	print_current_time("Second attempt to initiate the progressive assembly: trimming the reads ends, decreasing the k-mer size, and re-running the assembly");
	print "Second attempt to initiate the progressive assembly: trimming the reads ends, decreasing the k-mer size, and re-running the assembly\n";

      SECOND_ATTEMPT:
	foreach my $value (@decreasing_values) {
	    my $decreased_kmer = $kmer - $value;

	    open(NEW_FILE, ">$output/sequences_trimmed.fasta") or die "ERROR: Could not create file $output/sequences_trimmed.fasta : $!\n";
	    open(HANDLE, $input_file_for_assembler) or die "ERROR: Could not open file $input_file_for_assembler : $!\n";

	    while(my $line = <HANDLE>) {
		if ($line =~ m/^(>.*)\n/) {
		    print NEW_FILE $1, "\n";
		}
		else {
		    print NEW_FILE substr ( $line, 0, ( (length $line) -1) - $value);
		    print NEW_FILE "\n";
		}
	    }
	    close NEW_FILE;
	    close HANDLE;
	    if ($assembler_name =~ m/ABySS/) {
		system "abyss -k $decreased_kmer $assembler_parameters -o $output/contigs.fasta $output/sequences_trimmed.fasta >>$output/assembler.log 2>>$output/assembler.log";
		$stop_progressive_assembly = check_assembler_reconstruction("$output/contigs.fasta");
	    }
	    elsif ($assembler_name =~ m/Velvet/) {
		system "velveth $output/ $decreased_kmer $output/sequences_trimmed.fasta >> $output/assembler.log 2>> $log";
		system "velvetg $output/ $assembler_parameters >> $output/assembler.log 2>> $log";

		remove_intermediary_assembler_files();
		$stop_progressive_assembly = check_assembler_reconstruction("$output/contigs.fa");
	    }
	    elsif ($assembler_name =~ m/SOAPdenovo/) {
		my $soap_assembler;

		if ($kmer <= 127 && $kmer > 63) {
		    $soap_assembler = "SOAPdenovo-127mer";
		}
		elsif ($kmer <= 63 && $kmer > 31) {
		    $soap_assembler = "SOAPdenovo-63mer";
		}
		elsif ($kmer <= 31) {
		    $soap_assembler = "SOAPdenovo-31mer";
		}

		open(SOAP_CONF, ">$output/soap_configuration_file.temp") or die "ERROR: Could not create SOAP configuration file $output/soap_configuration_file.temp : $!\n";
		print SOAP_CONF "max_rd_len=1000000\n";
		print SOAP_CONF "[LIB]\n";
		print SOAP_CONF "asm_flags=1\n";
		print SOAP_CONF "f=$output/sequences_trimmed.fasta";
		close SOAP_CONF;

		system "$soap_assembler pregraph -K $decreased_kmer -s $output/soap_configuration_file.temp $assembler_parameters -o $output/contigs >> $output/assembler.log 2>> $log";
		system "$soap_assembler contig -g $output/contigs >> $output/assembler.log 2>>$log";

		remove_intermediary_assembler_files();
		$stop_progressive_assembly = check_assembler_reconstruction("$output/contigs.contig");
	    }
	    if(!$stop_progressive_assembly) {
		last SECOND_ATTEMPT;
	    }

	    #
	    # remove intermediate file
	    #
	    unlink "$output/sequences_trimmed.fasta";
	}
    }

    #
    # THIRD attempt
    #
    if ($assembler_name =~ m/(Velvet)|(SOAPdenovo)/i && $stop_progressive_assembly) {
	my $output_assembler;

	if ($assembler_name =~ m/Velvet/i) {
	    $output_assembler = "$output/contigs.fa";
	}
	else{
	    $output_assembler = "$output/contigs.contig";
	}
	print_current_time("Third attempt to initiate the progressive assembly: switching $assembler_name for ABySS");
	print "Third attempt to initiate the progressive assembly: switching $assembler_name for ABySS\n";

	system "abyss -k $kmer -c1 -e2 -o $output_assembler $input_file_for_assembler >>$output/assembler.log 2>>$output/assembler.log";
	$stop_progressive_assembly = check_assembler_reconstruction("$output/contigs.contig");

	if ($stop_progressive_assembly) {

	  THIRD_ATTEMPT:
	    foreach my $value (@decreasing_values) {
		my $decreased_kmer = $kmer - $value;

		system "abyss -k $decreased_kmer -c1 -e2 -o $output_assembler $input_file_for_assembler >>$output/assembler.log 2>>$output/assembler.log";
		$stop_progressive_assembly = check_assembler_reconstruction("$output_assembler");

		if(!$stop_progressive_assembly) {
		    last THIRD_ATTEMPT;
		}
	    }
	}
    }
    return($stop_progressive_assembly);
}

sub check_progress {
    my $last_contig_file = shift;
    my $next_contig_file = shift;
    my $seed_type = shift;
    my $output = shift;
    my $positive_contigs = shift;
    my $current_contigs = shift;
    my $round = shift || "final";

    #
    # Formats the subject and starts the similarity search 
    #
    my $stop_progressive_assembly = 0;

    if ($seed_type !~ m/Profile HMM/) {
	system "makeblastdb -in $next_contig_file -dbtype nucl -logfile $output/makeblastdb.log";

	if ($seed_type =~ m/DNA/) {
	    !system "blastn -query $last_contig_file -db $next_contig_file -perc_identity 90 -dust no -outfmt \"6 qseqid qlen length sseqid slen frames\" -out $output/check_progress.txt -max_hsps 1 2>>$log" or
		die "ERROR: Could not run blastn: $!\nCommand: blastn -query $last_contig_file -db $next_contig_file -perc_identity 90 -dust no -outfmt \"6 qseqid qlen length sseqid slen frames\" -out $output/check_progress.txt -max_hsps 1 2>>$log";
	}
	elsif ($seed_type =~ m/Protein/) {
	    !system "tblastn -query $last_contig_file -db $next_contig_file -seg no -outfmt \"6 qseqid qlen length sseqid slen frames\" -out $output/check_progress.txt -max_hsps 1 2>>$log" or
		die "ERROR: Could not run tblastn: $!\nCommand: tblastn -query $last_contig_file -db $next_contig_file -seg no -outfmt \"6 qseqid qlen length sseqid slen frames\" -out $output/check_progress.txt -max_hsps 1 2>>$log\n";
	}

	#
	# remove intermediate files
	#
	unlink glob "$next_contig_file.n*";
    }
    else {
	print "";
	system "transeq -sformat pearson $next_contig_file $next_contig_file.transeq -frame 6 2>>$log";
	system "hmmsearch $hmmsearch_parameters --tblout $output/check_progress.txt -o $output/check_progress_full.txt $last_contig_file $next_contig_file.transeq 2>>$log";

	#
	# remove intermediate files
	#
	unlink "$output/check_progress_full.txt";
	unlink "$next_contig_file.transeq";
    }

    open(RESULT, "<$output/check_progress.txt") or die "ERROR: Could not open RESULT file $output/check_progress.txt (2007): $!\n";

    #
    # maybe more than one frame of one contig can match the query sequence (it happens occasionally)
    # Since the DNA sequence is important for the progressive assembly process, we must avoid redundant contigs in GenSeed-HMM.
    #
    my %only_one_frame_per_contig;
    my %grew;

    while(my $line = <RESULT>) {
	chomp $line;

      	if ($seed_type !~ m/Profile HMM/) {
	    my @blast_result = split(/\t/, $line);

	    #
	    # to avoid saving the same contig more than once per round
	    #
	    if(exists $only_one_frame_per_contig{$blast_result[3]}) {
		next;
	    }
	    else {
		$only_one_frame_per_contig{$blast_result[3]}++;
	    }

	    #
	    # last contig (or original seed) file is incorporated in the next contig file (at least $align_threshold)?
	    # => (alignment length / last contig length) * 100
	    #
	    if ((($blast_result[2]/$blast_result[1])*100) >= $align_threshold) {
		$grew{$blast_result[0]}++;
		#
		# if new contig length is shorter than max_contig_length
		#
		if($blast_result[4] <= $max_contig_length) {
		    #
		    # if contig increase in size in comparison to the last round
		    #
		    if ($blast_result[4] > $blast_result[1]) {
			$positive_contigs->{"r$round.$blast_result[3]"} = return_sequence_from_fasta($next_contig_file, $blast_result[3], $blast_result[5]);
		    }
		    else {

			#
			# if seed-contig not extend extract the sequence, trim the ends, and run a local assembly process
			#
			print "Trimming and trying to extend contig $blast_result[3]... ";
			print_current_time("Trimming and trying to extend contig $blast_result[3]");

			my $sequence_to_be_trimmed = return_sequence_from_fasta($next_contig_file, $blast_result[3], $blast_result[5]);

			my $contig_length_to_be_overcome = $blast_result[4];

			#
			# the actual contigs can be smaller than the previous (WHY? Perhaps too many errors create a tangled de Bruijn graph 
			# and the algorithm performs an agressive trim which causes a reduction in contig size).
			#
			my $tries = 0;
			my $stop_trimming_process = 0;
			my $previous_contig_length = $blast_result[1];
			my $flag;

			while($tries < 3 && !$stop_trimming_process) {
			    my $trimmed_sequence = plain_trim($sequence_to_be_trimmed);
			    last if $trimmed_sequence eq "-1";

			    ($stop_trimming_process) = trimming_and_extending($positive_contigs,
									      $trimmed_sequence,
									      $contig_length_to_be_overcome,
									      $previous_contig_length,
									      \$flag,
									      $output,
									      $tries,
									      $round);
			    $sequence_to_be_trimmed = $trimmed_sequence;
			    $tries++;
			}

			if (($error_hash{'SOAPdenovo-31mer'} || $error_hash{'SOAPdenovo-63mer'} || $error_hash{'SOAPdenovo-127mer'}) && $assembler_name =~ m/SOAPdenovo/i) {
			    print "ERROR: Could not find SOAPdenovo in the command path. Select another assembler.\n";
			    exit;
			}
			#
			# if stop trimming process is zero, contig didn't extend; put the contig in the final contig file
			#
			print "Done\n";

			if (!$stop_trimming_process) {
			    open(FINAL_FILE, ">>$output/final_result_dir/final_contigs.fasta") or die "ERROR: Could not write to file $output/final_result_dir/final_contigs.fasta (2083): $!\n";
			    my $seq;

			    if (!$flag) {
				$seq = return_sequence_from_fasta($next_contig_file, $blast_result[3], $blast_result[5]);
			    }
			    else {
				$seq = return_sequence_from_fasta($last_contig_file, $blast_result[0], $blast_result[5]);
			    }

			    my $length = length $seq;

			    print "Contig final_contig_$round.$blast_result[3] reconstructed with $length bp - finishing its progressive assembly (trimming did not improve reconstruction)\n";
			    print_current_time("Contig final_contig_$round.$blast_result[3] reconstructed with $length bp - finishing its progressive assembly (trimming did not improve reconstruction))");

			    print FINAL_FILE ">final_contig_$round.$blast_result[3]\n";
			    print FINAL_FILE $seq, "\n";
			    close FINAL_FILE;

			    delete $current_contigs->{$blast_result[0]} unless $current_contigs eq 'last';
			}
		    }
		}
		else {
		    my $seq = return_sequence_from_fasta($next_contig_file, $blast_result[3], $blast_result[5]);

		    print "Contig r$round.$blast_result[3] reached the maximum number of bases $blast_result[4] bp - finishing its progressive assembly\n" unless $current_contigs eq 'last';
		    print_current_time("Contig r$round.$blast_result[3] reached the maximum number of bases $blast_result[4] bp - finishing its progressive assembly") unless $current_contigs eq 'last';

		    if($current_contigs eq 'last') {
			$positive_contigs->{"r$round.$blast_result[3]"} = $seq;
		    }
		    else {
			open(FINAL_FILE, ">>$output/final_result_dir/final_contigs.fasta") or die "ERROR: Could not write to file $output/final_result_dir/final_contigs.fasta (2107): $!\n";
			print FINAL_FILE ">r$round.$blast_result[3]\n";
			print FINAL_FILE $seq, "\n";
			close FINAL_FILE;
			delete $current_contigs->{$blast_result[0]};
		    }
		}
	    }
	    else {

		#
		# Due to limitations in Newbler, we need to split the contigs bigger than 2000 to run the assembly. But when we do this, we introduce a drawback
		# in the progressive assembly. If we have regions with repetitions the product of the assembly process can be fragmented and the last contig seed
		# will not be present in the former contigs. So when this happens we need to finish the progressive assembly process.
		#
		if (defined $round && looks_like_number($round) && $round == 1) {
		    print "Seed $blast_result[0] not found in reconstructed contig. Stopping progressive assembly for contig $blast_result[3].\n";
		    print_current_time("Seed $blast_result[0] not found in reconstructed contig. Stopping progressive assembly for contig $blast_result[3].");
		    next;
		}

		if ($current_contigs =~ m/last/) {
		    next;
		}
		unless(exists $grew{$blast_result[0]}) {
		    open(FINAL_FILE, ">>$output/final_result_dir/final_contigs.fasta") or die "ERROR: Could not write to file $output/final_result_dir/final_contigs.fasta (2133): $!\n";
		    my $seq = return_sequence_from_fasta($last_contig_file, $blast_result[0], $blast_result[5]);

		    my $length = length $seq;
		    my $r = $round - 1;

		    print "Contig $blast_result[0] reconstructed with $length bp - finishing its progressive assembly (contig did not grow); New contig: $blast_result[4] (alignment length: $blast_result[2])\n";
		    print_current_time("Contig $blast_result[0] reconstructed with $length bp - finishing its progressive assembly (contig did not grow); New contig: $blast_result[4] (alignment length: $blast_result[2])");

		    print FINAL_FILE ">$blast_result[0]\n";
		    print FINAL_FILE $seq, "\n";
		    close FINAL_FILE;
		}

		delete $current_contigs->{$blast_result[0]} unless $current_contigs eq 'last';

		if (($error_hash{'SOAPdenovo-31mer'} || $error_hash{'SOAPdenovo-63mer'} || $error_hash{'SOAPdenovo-127mer'}) && $assembler_name =~ m/SOAPdenovo/i) {
		    print "ERROR: Could not find SOAPdenovo in the command path. Select another assembler.\n";
		    exit;
		}
	    }
	    delete $current_contigs->{$blast_result[0]} unless $current_contigs eq 'last';
	}
	else {
	    if ($line !~ /^#/) {
		my @hmmsearch_result = split (/\s+/, $line);
		if ($hmmsearch_result[0] =~ m/(\S+)(_\d)/) {
		    my $tmp = $1;
		    if(!exists $only_one_frame_per_contig{$tmp}) {
			$positive_contigs->{"r$round.$tmp"} = return_sequence_from_fasta($next_contig_file, $tmp);
			$only_one_frame_per_contig{$tmp}++;
		    }
		}
	    }
	}
    }
    close RESULT;

    if (!%$positive_contigs) {
	$stop_progressive_assembly = 1;
	if ($current_contigs !~ m/last/) {
	    print_current_time("No contig has the necessary specifications to continue the assembly process");
	    print "No contig has the necessary specifications to continue the assembly process.\n";
	}
	else {
	    print "Done.\n";
	}
    }
    return ($stop_progressive_assembly);
}

sub return_sequence_from_fasta {
    my $fasta_filename = shift;
    my $fasta_header = shift;
    my $frames = shift;

    local $/ = undef;
    open(FASTA, $fasta_filename) or die "ERROR: Could not open file $fasta_filename (2406): $!\n";

    my $fasta_sequence;

    while(my $line = <FASTA>) {
	if ($line =~ m/(>)($fasta_header\n)((A|T|C|G|N|\n)+)/i || $line =~ m/(>)($fasta_header .*?\n)((A|T|C|G|N|\n)+)/i) {
	    $fasta_sequence = $3;
	}
	else {
	    return -1;
	}
    }
    $fasta_sequence =~ s/\n//g;

    #
    # verifies the next contig file orientation and if they're reverse complementar of each other
    # take this in consideration and correct the contig sequence for the next round
    #
    if(defined $frames) {
	my ($id_orientation, $subject_orientation) = split(/\//, $frames);
	if ($id_orientation != $subject_orientation) {
	    $fasta_sequence = reverse($fasta_sequence);
	    $fasta_sequence =~ tr/ATCGatcg/TAGCtagc/;
	}
    }
    return($fasta_sequence);
}

sub plain_trim {
    my $sequence_for_trim = shift;

    #
    # abort trimming if $sequence_for_trim is < 2.5 * $next_seed_size
    #
    if(length($sequence_for_trim) < 2.5 * $next_seed_size) {
	print " Contig could not be trimmed: shorter than 2.5 * ext_seed_size ($next_seed_size). ";
	print_current_time("Contig could not be trimmed: shorter than 2.5 * ext_seed_size ($next_seed_size)\n");
	return -1;
    }

    #
    # Trim the contig by 0.25 x the size of the extension seed
    #
    my $trimmed_sequence;
    my $trimming_value = 0.25;
    my $offset = int($next_seed_size * $trimming_value);

    $trimmed_sequence = substr($sequence_for_trim, $offset, -$offset);
    return($trimmed_sequence);
}

sub trimming_and_extending {
    my $positive_contigs = shift;
    my $trimmed_sequence = shift;
    my $contig_length_to_be_overcome = shift;
    my $previous_contig_length = shift;
    my $flag = shift;
    my $output = shift;
    my $tries = shift;
    my $round = shift;

    #
    # creating temporary contig file and temporary seed files for similarity search
    #
    my $next_seed_file = "$output/trimmed_seed.fasta";
    my $contig_file = "$output/trimmed_contig.fasta";
    my $blast_result = "$output/trimmed_blast.txt";
    my $recruited_reads_list = "$output/trimmed_read_list.txt";
    my $recruited_reads_filename = "$output/trimmed_read_list.fasta";
    my $new_blast_result = "$output/trimmed_newblast.txt";

    my $new_contig_file;
    my $new_contig_quality_file;
    my @recruited_reads;

    open(TEMPORARY_CONTIG, ">$contig_file") or die "ERROR: Could not create file $contig_file (2262): $!\n";

    print TEMPORARY_CONTIG ">trimmed_contig\n";
    print TEMPORARY_CONTIG $trimmed_sequence;

    close TEMPORARY_CONTIG;

    generate_sequence_for_similarity_search($contig_file, $next_seed_file);

    #
    # running similarity search, filtering and retrieving the blast hits
    #
    !system "blastn -query $next_seed_file -db $db $blastn_parameters -outfmt \"6 sseqid qlen length\" -out $blast_result 2>>$log" or
	die "ERROR: Could not run blastn : $!\nCommand: blastn -query $next_seed_file -db $db $blastn_parameters -outfmt \"6 sseqid qlen length\" -out $blast_result 2>>$log\n";


    open(BLAST_RESULT, $blast_result) or die "ERROR: Could not open file $blast_result : $!\n";

    while (my $line = <BLAST_RESULT>) {
	my @filter_blast = split(/\t/, $line);
	if ((($filter_blast[2]/$filter_blast[1])*100) >= $min_hsp_length) {
	    push (@recruited_reads, $filter_blast[0]);
	}
    }
    close BLAST_RESULT;

    open (FILE, ">$recruited_reads_list") or die "ERROR: Could not open file $recruited_reads_list : $!\n";

    foreach my $hit_name (@recruited_reads) {
        print FILE $hit_name;
        print FILE "\n";
    }

    close FILE;
    close BLAST_RESULT;

	my $cmd1 = "grep -F --no-group-separator -A 1 -f $recruited_reads_list $db > $recruited_reads_filename";
	`$cmd1`;

    #
    # fix recruited_reads sequence IDs (remove lcl| -- or whatever 3 letters appear after > and before |)
    #
    system("sed -i 's/>...\|/>/' $recruited_reads_filename");

    #
    # run the assembly process and remove intermediate assembly files
    #
    if ($assembler_name =~ m/(Velvet)/i) {
	run_velvet($contig_file, $recruited_reads_filename, \$new_contig_file, $round);
    }
    elsif ($assembler_name =~ m/SOAPdenovo/i) {
	run_soapdenovo($contig_file, $recruited_reads_filename, \$new_contig_file, $round);
    }
    elsif ($assembler_name =~ m/ABySS/i) {
	run_abyss($contig_file, $recruited_reads_filename, \$new_contig_file, $round);
    }
    elsif ($assembler_name =~ m/Newbler/i) {
	run_newbler($contig_file, $recruited_reads_filename, \$new_contig_file, \$new_contig_quality_file, $round);
	unlink "$output/454AllContigs.qual" if (-e "$output/454AllContigs.qual");
    }
    else {
	run_cap3($contig_file, $recruited_reads_filename, \$new_contig_file, \$new_contig_quality_file);
	system "rm $output/input_sequences.fasta.cap.contigs.qual" if (-e "$output/input_sequences.fasta.cap.contigs.qual");
    }

    #
    # identify positive contigs in this assembly, check if the former contig is bigger than the original, and remove intermediate files
    #
    unlink "$next_seed_file";
    unlink "$blast_result";
    unlink "$recruited_reads_list";
    unlink glob "$recruited_reads_filename*";

    if ($new_contig_file) {
	system "makeblastdb -in $new_contig_file -dbtype nucl -logfile $output/makeblastdb.log";
	system "blastn -query $contig_file -max_target_seqs 1 -max_hsps 1 -db $new_contig_file -outfmt \"6 qseqid qlen length sseqid slen frames\" -out $new_blast_result 2>>$log";
	unlink "$contig_file";

	open(RESULT_BLAST, $new_blast_result);

	my $rand = int(10000 * rand());

	while (my $line = <RESULT_BLAST>) {
	    my @blast_result = split(/\t/, $line);
	    if ((($blast_result[2]/$blast_result[1])*100) >= 95) {
		if ($blast_result[4] > $contig_length_to_be_overcome) {
		    if ($blast_result[4] > $previous_contig_length) {
			$positive_contigs->{"r$round.te.$blast_result[3].$rand"} = return_sequence_from_fasta($new_contig_file, $blast_result[3], $blast_result[5]);
			unlink glob "$new_contig_file*";
			unlink "$new_blast_result";
			close RESULT_BLAST;
			$$flag = 0;
			return(1);
		    }
		    else {

			#
			# the actual contigs can be smaller than previous so in the final contig files will print the bigger contig. 
			# This flag does this, telling the program how contig will go to the final contig file
			#
			$$flag = 1;
		    }
		}
	    }
	    else {
		$$flag = 1;
	    }
	}
	unlink glob "$new_contig_file*";
	unlink "$new_blast_result";
	close RESULT_BLAST;
	return(0);
    }
    else {
	unlink "$contig_file";
	return(0);
    }
}

sub get_file {
    my $recruited = shift;
    my $last_round = shift;
    my $assembler = shift;
    my $next_round_qual = shift;
    my $round = shift;

    my $input;

    if($seed_information->{'type'} =~ m/DNA/) {
       	$input = "$output/input_sequences.$ext";
        if($assembler eq "newbler") {

	    #
	    # If no SFF files were given, put recruited reads (FASTA format) in input file for assembler, or two files (FASTA and FASTA + Qual);
	    # otherwise, $recruited_reads will be a new SFF file with the reads (flowgram format)
	    #
	    if($db_info->{'db_type'} !~ m/sff/i) {
		system "cat $recruited > $input";

		if(defined $db_info->{'retr_qual'}) {
		    system "cp $db_info->{'retr_qual'} $input\.qual";
		}
	    }
	    if($round > 1 || $db_info->{'db_type'} =~ m/sff/i) { 
		generate_newbler_input_file($last_round, $input, $next_round_qual);
	    }
	}
	else {
	    system "cat $recruited $last_round > $input";
	}
	return($input);
    }
    else {
	return($recruited);
    }
}

sub fastq2fastaseq {
    my $fq_file = shift;

    my $fa_file = $fq_file . ".fasta";

    my $fa_qual_file = $fq_file . ".fasta.qual";

    print "\nConverting FASTQ file $fq_file to FASTA and QUAL... \n";
    print_current_time("Converting FASTQ file $fq_file to FASTA and QUAL... ");
    
    open(FQ, "<$fq_file") 
	or die "ERROR: Could not open FASTQ file $fq_file (2457): $!\n";
    open(FA, ">$fa_file") 
	or die "ERROR: Could not create FASTA file $fa_file (2459): $!\n";
    open(FAQ, ">$fa_qual_file") 
	or die "ERROR: Could not create QUALITY file $fa_qual_file (2461): $!\n";

    my $count = 0;

    while(<FQ>) {
	if($count % 4 == 0) {
	    my $display = $count / 4;
	    print "\r\tConverted $display reads..." if $display % 10000 == 0 || $display < 10000;
	    print FA ">", substr($_, 1);
	    print FAQ ">", substr($_, 1);
	}
	elsif($count % 4 == 1) {
	    print FA;
	}
	elsif($count % 4 == 3) {
	    chomp;
	    my @Qualities = split('', $_);
	    foreach my $quality (@Qualities) {
		print FAQ ord($quality) - $offset.' ';
	    }
	    print FAQ "\n";
	} 
	$count++;
    }
    my $display = $count / 4;
    print "\rDone ($display reads).                                                   \n";
    close FAQ;
    close FA;
    return $fa_file;
}

sub quality_splitter {
    my $file = shift;

    my $sfile = $file . "_split";
    my $size = 1500;
    my $overlap = 1300;
    my $flank = $size - $overlap;

    open(QF, "<$file") or die "ERROR: Could not open quality file $file (2783): $!\n";
    open(SQF, ">$sfile") or die "ERROR: Could not create split quality file $sfile (2784): $!\n";
    local $/ = '>';

    while(<QF>) {
	chomp;
	my($id, $quals);
	if($_ =~ /(\S+).*?\n(.+)/s) {
	    $id = $1;
	    $quals = $2;
	}
	else {
	    next;
	}

	my(@quals) = split(/\s+/, $quals);
	my $start = 0;
	my $length = scalar(@quals) - 1;
	my $end = 0;

	for(my $start = 0; $end < $length; $start += $flank) {
	    $end = ($start + $size - 1) > $length ? $length : ($start + $size - 1);
	    my $pstart = $start + 1;
	    my $pend = $end + 1;
	    my $data = join(" ", @quals[$start..$end]);
	    print SQF ">$id\_$pstart\-$pend\n$data\n";
	}
	
    }
    close SQF;
    return $sfile;
}

sub FA_index {
    my $input = shift;
    my $table;

    print "Indexing file $input...\n";
    print_current_time("Indexing file $input...");

    my $fh = create_FH();
    my $ifh = create_FH();
    open($ifh, "<$input") or die "ERROR: Could not open file $input for indexing (2799): $!\n";
    my $base = $input . ".idx";
    my $base2 = $input . ".iddb";

    if(-e $base && -e $base2) {
	print "Input file $input has already been indexed (files $base and $base2 already present).\n";
	print_current_time("Input file $input has already been indexed (files $base and $base2 already present).");

	open($fh, "+<$base") or die "ERROR: Could not open index file $base (2810): $!\n";
	my $st = time();
	$table = retrieve($base2);
	my $end = time();
	my $diff = $end - $st;
	print "Retrieval of ID table from disk performed in $diff seconds.\n";
	print_current_time("Retrieval of ID table from disk performed in $diff seconds.\n");
	return $fh, $ifh, $table;
    }
    else {
	unlink $base if(-e $base);
	unlink $base2 if(-e $base2);

	open($fh, "+>$base") or die "ERROR: Could not create index file $base (2826): $!\n";
	my $cnt = 0;

	local $/ = ">";
	my $re = qr/(\S+).*?\n.+/s;
	<$ifh>;
	my $offset = tell($ifh);
	my $st = time();

	while(<$ifh>) {
	    my $Fid;
	    if($_ =~ /$re/o) {
		$Fid = $1;
	    }
	    die "ERROR: duplicate sequence identifier $Fid\n" if exists $$table{$Fid};
	    $$table{$Fid} = $cnt;
	    my $old_off = $offset;
	    $offset = tell($ifh);
	    $cnt++;
	    print $fh pack("QL", $old_off, $offset - $old_off - 1);
	}
	my $end = time();
	my $diff = $end - $st;
	print "Indexing of $cnt entries performed in $diff seconds.\n";
	print_current_time("Indexing of $cnt entries performed in $diff seconds.");

	print "Storing of ID table to disk...";
	$st = time();
	store $table, $base2;
	$end = time();
	$diff = $end - $st;
	print " performed in $diff seconds.\n";
	print_current_time("Storing of ID table to disk performed in $diff seconds.\n");
    }
    return $fh, $ifh, $table;
}

sub get_indexed_recs {
    my $fh = shift;
    my $ifh = shift;
    my $query = shift;
    my $outfile = shift;
    my $table = shift;

    open(Q, "<$query") or die "ERROR: Could not open query file $query (2866): $!\n";
    open(ORR, ">$outfile") or die "ERROR: Could not create recruited reads file $output (2859): $!\n";
    my $size = length(pack("QL", 0, 0));

    while(my $FID = <Q>) {
	chomp $FID;
	my($content, $data);
	my $index_offset = $$table{$FID} * $size;
	seek($fh, $index_offset, 0) or die "ERROR: Could not seek to index position $index_offset (2872): $!\n";
	read($fh, $content, $size);
	my($data_offset, $datalen) = unpack("QL", $content);
	seek($ifh, $data_offset, 0) or die "ERROR: Could not seek to data position $data_offset (2875): $!\n";
	read($ifh, $data, $datalen);
	print ORR ">$data";
    }
    print ORR "\n";
    close ORR;
}

sub create_FH {
    local *FH;
    return *FH;
}
