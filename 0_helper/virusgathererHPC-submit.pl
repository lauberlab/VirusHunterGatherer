#!/usr/bin/perl

use warnings;
use strict;
#use Data::Dumper;

use lib "/home/lauber/lib/perl";
use Parallel::ForkManager;

use lib "/home/lauber/lauber-2015-virusHunter";
use      virusHGutils;


# parameters
my $memPerCPUpp     = 1875;
my $memPerCPUass    = 2583;
my $cpusPerTask     = 12;
my @ncbiIPs         = ("130.14.250.7","130.14.250.10","130.14.250.11",
                       "130.14.250.12","130.14.250.13");
my @ebiIPs          = ("fasp.sra.ebi.ac.uk");
my $useNCBI         = 1;
my $useEBI          = 1;
my $doClean         = 1;
my $combineByTxid   = 0;
my $usrAssemblyName = "";
my $useBlastSeed    = 0;
my $userSeed        = "";
my $useHitReads     = "";
my $copyResults     = 1;
my $readQcut        = 15;
my $ppOnly          = 0;
my $GS_assembler    = "newbler";
my $GS_seed_size    = 45;
my $GS_blast_E      = 0.001;
my $GS_blast_ide    = 85;
my $GS_blast_ws     = 11; # corresponds to blastn task: 7 (blastn-short), 11 (blastn), 28 (megablast)
my $GS_contig_len   = 10000;
my $GS_max_itr      = 500;
my $GS_aln_th       = 25;
my $GS_use_qual     = "no";
my $GS_starter      = 0;
my $GS_cap_identity = 85;
my $GS_cap_overlap  = 20;
my $GS_cap_overhang = 99;
my $GS_direction    = "b";
my $maxHoursPP      = 24;
my $maxHoursAss     = 48;
my $TNT             = ""; # assemble full data set using Trinity or SPAdes instead of Genseed
my $useLustre       = 0;
my $userOutdir      = "";
my $wait            = 1;
my $fastq_input     = 0;
my $userURL         = 0;
my $use_reservation = "";

# input
if ( $#ARGV < 3 ){
die("\n # # # # # # # # # # # # # # # # # # # # # # # # #
 # virusgatherer HPC version 2.9                 #
 #  assembly of viral sequences from SRA runs    #
 #  using GenSeed-HMM as core engine             #
 #  or using Trinity / SPAdes                    #
 # designer:  Stefan Seitz and Chris Lauber      #
 # developer: Chris Lauber                       #
 # email:     chris.lauber\@tu-dresden.de         #
 # # # # # # # # # # # # # # # # # # # # # # # # #\n
usage: virusgathererHPC-submit.pl <family> <sra_accs_file> <project> <tasks> [options]\n
options:
\t-t=<integer>\tuse <t> threads per task (default: $cpusPerTask)
\t-m\t\tmerge SRA runs with same NCBI taxon ID (default: each SRA run separately)
\t-x\t\tdo not delete temporary files
\t-b\t\tuse blastp/n for initial seed search (default: hmmsearch)
\t-nowait\t\tdo not wait for completion of all assemblies
\t-v=<project>\tuse contigs produced by virushunter search as seed
\t\t\t(default: initial seed search by hmmsearch)
\t-se=<string>\tfile with user-provided seed sequence(s)
\t-c\t\tdo not copy assembly results to project folder (default: copy)
\t-c=<string>\tcombine all SRA runs into assembly <c>
\t-o=<string>\tsave result in output directory <o>
\t-ha=<integer>\tmaximum runtime of the assembly step in hours (default: $maxHoursAss)
\t-q=<integer>\tquality cut-off for read trimming (default: $readQcut)
\t-a=<string>\tGenSeed-HMM assembler parameter (default: $GS_assembler)
\t-u\t\tuse auxiliary starter assembler (cap3) for newbler (default: none)
\t-s=<integer>\tGenSeed-HMM seed size parameter (default: $GS_seed_size)
\t-e=<float>\tGenSeed-HMM blastn E-value parameter (default: $GS_blast_E)
\t-i=<integer>\tGenSeed-HMM blastn identity parameter (default: $GS_blast_ide)
\t-w=<integer>\tGenSeed-HMM blastn word size parameter (default: $GS_blast_ws)
\t-l=<integer>\tGenSeed-HMM max contig length parameter (default: $GS_contig_len)
\t-k=<integer>\tGenSeed-HMM max iterations parameters (default: $GS_max_itr)
\t-h=<integer>\tGenSeed-HMM alignment threshold parameter (default: $GS_aln_th)
\t-p=<integer>\tGenSeed-HMM CAP3 percent identiy parameter (default: $GS_cap_identity)
\t-ol=<integer>\tGenSeed-HMM CAP3 absolute overlap parameter (default: $GS_cap_overlap)
\t-oh=<integer>\tGenSeed-HMM CAP3 percent overhang parameter (default: $GS_cap_overhang)
\t-di=<l|r|b>\tGenSeed-HMM expansion direction (default: $GS_direction)
".
#\t-q=<yes|no>\tGenSeed-HMM use quality parameters (default: $GS_use_qual)
"\t-r=<integer>\tuse <r> Mb of memory for assembly job (default: $memPerCPUass)
\t\t\teither of [ 1875 | 2583 | 3875 | 5250 | 7875 | 10583 | 36500 ]
\t-ppOnly\t\tonly download and preprocess the data (default: subsequent assembly)
\t-ncbi\t\tdownload only from NCBI (default: from EBI and NCBI)
\t-ebi\t\tdownload only from EBI  (default: from EBI and NCBI)
\t-TNT=<s|t>\tassemble full data set using SPAdes (s, under construction) or
\t\t\tTrinity(t, only recommended for RNA-seq data)
\t-lustre\t\tsave files on /lustre/ssd instead of /scratch (10 SRAs at most)
\t-fastq\t\tprovided are not SRA identifiers for download, but directory names
\t\t\tof local fastq.gz files
\t-url\t\tuse path from RunInfo for download (deprecated)
\t-rsv=<string>\tuse reservation <rsv> on Taurus (default: submit to normal queue)
\n");
}

my $family      = shift;
my $sraFile     = shift;
my $project     = shift;
my $cores       = shift;
while ( $#ARGV >= 0 ){
	my $option = shift;
	if ( $option =~ /-t=(\d+)/ ){  $cpusPerTask     = $1; }
	if ( $option eq "-m"       ){  $combineByTxid   =  1; }
	if ( $option eq "-x"       ){  $doClean         =  0; }
	if ( $option eq "-b"       ){  $useBlastSeed    =  1; }
	if ( $option eq "-nowait"  ){  $wait            =  0; }
	if ( $option =~ /-v=(.+)/  ){  $useHitReads     = $1; }
	if ( $option =~ /-se=(.+)/ ){  $userSeed        = $1; }
	if ( $option =~ /-ha=(\d+)/){  $maxHoursAss     = $1; }
	if ( $option eq "-c"       ){  $copyResults     =  0; }
	if ( $option =~ /-c=(.+)/  ){  $usrAssemblyName = $1; }
        if ( $option =~ /-q=(\d+)/ ){  $readQcut        = $1; }
        if ( $option =~ /-a=(.+)/  ){  $GS_assembler    = $1; }
	if ( $option eq "-u"       ){  $GS_starter      =  1; }
        if ( $option =~ /-s=(\d+)/ ){  $GS_seed_size    = $1; }
        if ( $option =~ /-e=(.+)/  ){  $GS_blast_E      = $1; }
        if ( $option =~ /-i=(\d+)/ ){  $GS_blast_ide    = $1; }
        if ( $option =~ /-w=(\d+)/ ){  $GS_blast_ws     = $1; }
        if ( $option =~ /-l=(\d+)/ ){  $GS_contig_len   = $1; }
        if ( $option =~ /-k=(\d+)/ ){  $GS_max_itr      = $1; }
        if ( $option =~ /-h=(\d+)/ ){  $GS_aln_th       = $1; }
        if ( $option =~ /-p=(\d+)/ ){  $GS_cap_identity = $1; }
        if ( $option =~ /-ol=(\d+)/){  $GS_cap_overlap  = $1; }
        if ( $option =~ /-oh=(\d+)/){  $GS_cap_overhang = $1; }
        if ( $option =~ /-di=(.+)/){   $GS_direction    = $1; }
#        if ( $option =~ /-q=(.+)/  ){  $GS_use_qual     = $1; }
        if ( $option =~ /-r=(\d+)/ ){  $memPerCPUass    = $1; $memPerCPUpp = $memPerCPUass; }
        if ( $option eq "-ppOnly"  ){  $ppOnly          =  1; }
        if ( $option eq "-ncbi"    ){  $useEBI          =  0; }
        if ( $option eq "-ebi"     ){  $useNCBI         =  0; }
        if ( $option =~ /-TNT=(.+)/){  $TNT             = $1; }
        if ( $option eq "-lustre"  ){  $useLustre       =  1; }
	if ( $option =~ /-o=(.+)/  ){  $userOutdir      = $1; }
	if ( $option eq "-fastq"   ){  $fastq_input     =  1; }
	if ( $option eq "-url"     ){  $userURL         =  1; }
        if ( $option =~ /-rsv=(.+)/){  $use_reservation = $1; }
}

# global parameters
#my $bindir  = "/home/lauber/lauber-2015-virusHunter/developer";
my $bindir     = "/home/lauber/lauber-2015-virusHunter";
my $tmpdir0    = sprintf "%s/lauber-2015-virusGatherer", virusHGutils::get_workspace_path( "scratch" );
   $tmpdir0    = sprintf "%s/virusgatherer",             virusHGutils::get_workspace_path( "lustre"  )  if $useLustre == 1;
my $tmpdir     = "$tmpdir0/$family/$project";
   $tmpdir     = $userOutdir  if $userOutdir ne "";
my $resdir     = "/projects/p_sra/$family/results";
my %jobs       = ();

# create temporary family directory
if ( ! -d "$tmpdir0/$family" ){
	mkdir( "$tmpdir0/$family" );
}

# create temporary directory
if ( ! -d "$tmpdir" ){
	mkdir( "$tmpdir" );
}

# create project dir if not already present
if ( ! -d "$resdir/$project"  and  $copyResults != 0 ){
	mkdir( "$resdir/$project" );
}

# create assembly result dir if not already present
if ( ! -d "$resdir/$project/assemblies"  and  $copyResults != 0 ){
	mkdir( "$resdir/$project/assemblies" );
}

# read SRA IDS from file and determine
# taxon and NCIB txid
if ( $fastq_input ){
	print "[virusgatherer] reading fastq file info\n";
}else{
	print "[virusgatherer] getting taxonomy and run information from NCBI/SRA\n";
}
my $sras = {};
my $assemblies = {};
my ($ipiNCBI, $ipiEBI, $si)  = ( -1, -1, 0 );
my $provider = "";
open( SRA,"<$sraFile");
while ( my $line = <SRA> ){
	chomp($line);  next if $line eq "";
	$sras->{$line} = {};
	# provider for data download
	if ( $provider eq "" ){
		if ( $useNCBI == 1 ){  $provider = "ncbi"; }
		else{                  $provider = "ebi";  }
	}else{
		if    ( $provider eq "ncbi"  and  $useEBI  == 1 ){  $provider = "ebi";  }
		elsif ( $provider eq "ebi"   and  $useNCBI == 1 ){  $provider = "ncbi"; }
	}
	$sras->{$line}->{'provider'} = $provider;
	# IP/hostname for data download
	if ( $provider eq "ncbi" ){
		$ipiNCBI++;
		$ipiNCBI = 0 if $ipiNCBI > $#ncbiIPs;
		$sras->{$line}->{'ip'} = $ncbiIPs[ $ipiNCBI ];
	}else{
		$ipiEBI++;
		$ipiEBI = 0 if $ipiEBI > $#ebiIPs;
		$sras->{$line}->{'ip'} = $ebiIPs[ $ipiEBI ];
	}
	$si++;
	$sras->{$line}->{'ni'} = $si;
	# get taxon info from NCBI SRA
	$sras->{$line}->{'taxon'}  = "";
	$sras->{$line}->{'txid'}   = "";
	$sras->{$line}->{'sizeMB'} = "";
	if ( $fastq_input == 0 ){
		my $wget_cmd = "wget -q -O - 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=$line'";
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
				$sras->{$line}->{'txid'}   = $val[ $wget_i ] if ( $hdr[ $wget_i] eq "TaxID" );
				$sras->{$line}->{'sizeMB'} = $val[ $wget_i ] if ( $hdr[ $wget_i] eq "size_MB" );
				$sras->{$line}->{'URL'}    = $val[ $wget_i ] if ( $hdr[ $wget_i] eq "download_path" );
				$sras->{$line}->{'taxon'}  = $val[ $wget_i ] if ( $hdr[ $wget_i] eq "ScientificName" );
				$sras->{$line}->{'taxon'} =~ s/ /_/g;
				$sras->{$line}->{'taxon'} =~ s/\//_/g;
				$sras->{$line}->{'taxon'} =~ s/\(//g;
				$sras->{$line}->{'taxon'} =~ s/\)//g;
				$sras->{$line}->{'taxon'} =~ s/\'//g;
			}
		}
	}
	# sort SRA runs by assembly
	my $ass_name = $line;
	if ( $sras->{$line}->{'txid'} ne ""  and  $combineByTxid ){
		$ass_name  = "txid".$sras->{$line}->{'txid'};
	}
	if ( $sras->{$line}->{'taxon'} ne "" ){
		$ass_name .= "_".$sras->{$line}->{'taxon'};
	}
	if ( $fastq_input ){
		$line =~ /.*\/([^\/]+)$/;
		$ass_name = $1;
	}
	if ( $usrAssemblyName ne "" ){
		$ass_name  = $usrAssemblyName;
	}
	$sras->{$line}->{'assembly'} = $ass_name;
	$assemblies->{$ass_name} = [] if ! exists $assemblies->{$ass_name};
	push( @{ $assemblies->{$ass_name} }, $line );
}
close(SRA);

die( "\nRunning on /lustre/ssd/ only allowed for not more than 10 SRAs.\n\n" )  if ( $si > 10 and $useLustre == 1 );  
printf "[virusgatherer] will analyze $si SRA runs\n";
#print Dumper( $assemblies );


# fork into <cores> child processes and download data in parallel
if ( $fastq_input == 0 ){

print "[virusgatherer] downloading data from NCBI/SRA\n";
# init fork
my $pm = new Parallel::ForkManager($cores);
foreach my $sid ( sort{ $sras->{$a}->{'ni'} <=> $sras->{$b}->{'ni'} } keys %{$sras} ){
	# start fork
	my $pid = $pm->start and next;
	# download data
	my $ass = $sras->{$sid}->{'assembly'};
	if     ( -e "$tmpdir/$ass/data/data.fasta"                      and $TNT eq "" ){
		printf "[virusgatherer] \tdata already present for '$ass'. resuming\n";
	}elsif ( -e "$tmpdir/$ass/data/$sid"."_1.trim.q$readQcut.fastq" and $TNT ne "" ){
		printf "[virusgatherer] \tdata already present for '$ass'. resuming\n";
	}
	else{
		my $ip   = $sras->{$sid}->{'ip'};
		my $prov = $sras->{$sid}->{'provider'};
		my $cmd  = sprintf "/home/lauber/lauber-2015-virusHunter/virusgathererHPC-1-data.pl %s %s %s %s %s", $family, $project, $sid, $ip, $prov;
		   $cmd .= " -lustre"  if $useLustre == 1;
		   $cmd .= " -o=$userOutdir"  if $userOutdir ne "";
		   $cmd .= sprintf " -url=%s", $sras->{$sid}->{'URL'}  if $userURL == 1;
		`$cmd > $tmpdir/$sid-data.out 2> $tmpdir/$sid-data.err`;
	}
	# end fork
	$pm->finish;
}
# wait for all forks
$pm->wait_all_children;

}

# create assembly subdirectories and merge SRA runs if requested
print "[virusgatherer] initializing assembly directories and merging SRA runs if requested\n";
foreach my $ass ( keys %{ $assemblies } ){
	my $sid = $assemblies->{$ass}->[0];
	if     ( -e "$tmpdir/$ass/data/data.fasta"                      and $TNT eq "" ){
		next;
	}elsif ( -e "$tmpdir/$ass/data/$sid"."_1.trim.q$readQcut.fastq" and $TNT ne "" ){
		next;
	}
	else{
		mkdir( "$tmpdir/$ass" );
		mkdir( "$tmpdir/$ass/data" );
		foreach my $sraid ( @{ $assemblies->{$ass} } ){
			if ( $fastq_input ){
				`cp -r $sraid $tmpdir/$ass/data/`;
			}else{
				`mv $tmpdir/$sraid-*   $tmpdir/$ass/`;
				`mv $tmpdir/$sraid.sra $tmpdir/$ass/data/`;
			}
		}
	}
}


# submit assemblies
printf "[virusgatherer] submitting %d assembly jobs\n", scalar keys %{ $assemblies };
my $jobidfile = "$tmpdir/virusgathererHPC-$project-jobIDs.txt";
#open( JID, ">$jobidfile") or die("$0: could not open file: $!\n");
open( JID, ">>$jobidfile") or die("$0: could not open file: $!\n");
foreach my $ass ( keys %{ $assemblies } ){
	my $deps = "";
	my $sid  = $assemblies->{$ass}->[0];
	if     ( -e "$tmpdir/$ass/data/data.fasta"                      and $TNT eq "" ){
		;
	}elsif ( -e "$tmpdir/$ass/data/$sid"."_1.trim.q$readQcut.fastq" and $TNT ne "" ){
		;
	}
	else{
		foreach my $sraid ( @{ $assemblies->{$ass} } ){
			#my $ppJobId = submitPreprocessJob( $ass, $sraid, $readQcut, $cpusPerTask );
			if ( $fastq_input ){
				$sraid =~ /.*\/([^\/]+)$/;
				$sraid = $1;
			}
			my $ppJobId = submitPreprocessJob( $ass, $sraid, $readQcut, 1 ); # use only one CPU to avoid bug with autoadapt while splitting files
			print JID "$ppJobId\n";
			$deps .= ":".$ppJobId;
		}
	}
	if ( $ppOnly == 0 ){
		my $assJobId = submitAssemblyJob( $ass, $deps );
        	print JID "$assJobId\n";
	}
}
close(JID);


# wait for all jobs to finish
if ( $wait ){
	print "[virusgatherer] waiting for all jobs to finish\n";
	my %jobids = ();
	open( JIDS, "<$jobidfile" ) or die("$0: could not open file: $!\n");
	while (my $line = <JIDS> ){
		chomp($line);  next if $line eq "";
		$jobids{$line} = 1;
	}
	close(JIDS) or die("$0: could not close file: $!\n");
	while (1){
		last if ( scalar keys %jobids == 0 );
		foreach ( keys %jobids ){
			my $squeue = `squeue -u lauber | grep $_`;
			delete $jobids{$_} if $squeue eq "";
		}
		sleep(10);
	}

	# copy results to result directory
	if ( $copyResults ){
		foreach my $ass ( keys %{ $assemblies } ){
			if ( ! -d "$resdir/$project/assemblies/$ass" ){
				mkdir( "$resdir/$project/assemblies/$ass" );
				mkdir( "$resdir/$project/assemblies/$ass/results-$GS_assembler" );
				mkdir( "$resdir/$project/assemblies/$ass/log" );
			}
			`cp -r $tmpdir/$ass/results-$GS_assembler/* $resdir/$project/assemblies/$ass/results-$GS_assembler/`;
			`cp -r $tmpdir/$ass/*.out                   $resdir/$project/assemblies/$ass/log/`;
			`cp -r $tmpdir/$ass/*.err                   $resdir/$project/assemblies/$ass/log/`;
			`cp -r $tmpdir/$ass/*.txt                   $resdir/$project/assemblies/$ass/log/`;
			`cp -r $tmpdir/$ass/*.log                   $resdir/$project/assemblies/$ass/log/`;
		}
	}

	# delete temporary files
	if ( $doClean ){
		print "[virusgatherer] deleting temporary data and cleaning up\n";
		`rm -rf $tmpdir`;	
	}else{
		print "[virusgatherer] keeping temporary data at $tmpdir\n";
	}
}else{
	print "[virusgatherer] will not wait for jobs to complete\n";
}

# DONE
print "[virusgatherer] done\n";


# function to submit job for preprocessing an SRA run
sub submitPreprocessJob {
	my $assid   = shift;
	my $sraid   = shift;
	my $qcut    = shift;
	my $threads = shift;
	# which machines to use
	my $machines = "haswell,sandy";
	   $machines = "haswell" if ( $cpusPerTask > 12 );
	   $machines = "haswell" if ( $use_reservation ne "" );
	   $machines = "smp" if ( $memPerCPUass > 10583 );
	# compile job file for the search
	my $slurm = "#!/bin/bash\n";
	$slurm .= sprintf "#SBATCH --output=%s\n",		"$tmpdir/$assid/$sraid-preprocess.out";
	$slurm .= sprintf "#SBATCH --error=%s\n",		"$tmpdir/$assid/$sraid-preprocess.err";
	$slurm .= sprintf "#SBATCH --time=%d:00:00\n",		$maxHoursPP; # in hours
	$slurm .= sprintf "#SBATCH --nodes=%d\n",		1;
	$slurm .= sprintf "#SBATCH --ntasks-per-node=%d\n",	1;
	$slurm .= sprintf "#SBATCH --cpus-per-task=%d\n",	$threads;
	$slurm .= sprintf "#SBATCH --mem-per-cpu=%d\n",		$memPerCPUpp;
	$slurm .= sprintf "#SBATCH -p %s\n",			$machines;
	if ( $use_reservation ne "" ){
	$slurm .= sprintf "#SBATCH --reservation=%s\n",		$use_reservation;
	}
	$slurm .= sprintf "#SBATCH --job-name=%s\n\n",		"g$sraid-preprocess";
	$slurm .= sprintf "cd %s\n", $bindir;
	# load dependencies
	$slurm .= sprintf "module load Java/1.8.0_162\n";
	# parameters
	my $sraFile  = "$tmpdir/$assid/data/$sraid";
	   $sraFile .= ".sra"  if ( $fastq_input == 0 );
	$slurm .= sprintf "%s/virusgathererHPC-2-preprocess.pl %s -q=%d -t=%d",
                          $bindir, $sraFile, $qcut, $threads;
	$slurm .= " -tnt"     if ( $TNT ne "" );
	$slurm .= " -lustre"  if ( $useLustre   == 1 );
	$slurm .= " -fastq"   if ( $fastq_input == 1 );
	# create job file
	my $jobfile = "$tmpdir/$assid/$sraid-preprocess.slurm";

	printf "%s\n", $jobfile;
	
	open (SLURM,">$jobfile");
	print SLURM $slurm."\n";
	close(SLURM);
	# submit command
	my $cmd = "sbatch $jobfile | awk '{print \$4}'";
	# submit job file and log job ID
	my $jid = `$cmd`;  chomp($jid);
	# delete job file
	#`rm -rf $jobfile`;
	# return job ID
	return( $jid );
}


# function to submit job for an assembly project of one or several SRA runs
sub submitAssemblyJob {
	my $assid   = shift;
	my $depjobs = shift;
	# which machines to use
	my $machines = "haswell,sandy";
	   $machines = "haswell" if ( $cpusPerTask > 12 );
	   $machines = "haswell" if ( $use_reservation ne "" );
	   $machines = "smp" if ( $memPerCPUass > 10583 );
	#  compile job file for the search
	my $slurm = "#!/bin/bash\n";
	$slurm .= sprintf "#SBATCH --output=%s\n",		"$tmpdir/$assid/$assid-assembly.out";
	$slurm .= sprintf "#SBATCH --error=%s\n",		"$tmpdir/$assid/$assid-assembly.err";
	$slurm .= sprintf "#SBATCH --time=%d:00:00\n",		$maxHoursAss; # in hours
	$slurm .= sprintf "#SBATCH --nodes=%d\n",		1;
	$slurm .= sprintf "#SBATCH --ntasks-per-node=%d\n",	1;
	$slurm .= sprintf "#SBATCH --cpus-per-task=%d\n",	$cpusPerTask;
	$slurm .= sprintf "#SBATCH --mem-per-cpu=%d\n",		$memPerCPUass;
	$slurm .= sprintf "#SBATCH -p %s\n",			$machines;
	if ( $use_reservation ne "" ){
	$slurm .= sprintf "#SBATCH --reservation=%s\n",		$use_reservation;
	}
	$slurm .= sprintf "#SBATCH --job-name=%s\n\n",		"g$assid-assembly";
	$slurm .= sprintf "cd %s\n", $bindir;
	# load dependencies, if necessary
	#$slurm .= sprintf "module load python/3.6-anaconda4.4.0\n"  if $TNT eq "t";
	#$slurm .= sprintf "module load Python/3.6.4-intel-2018a\n"  if $TNT eq "t";
	$slurm .= sprintf "module load EMBOSS\n"                    if $TNT ne "t";
	#$slurm .= sprintf "module load Java/1.8.0_172\n"            if $TNT eq "t";
	$slurm .= sprintf "module load modenv/classic\n"            if $TNT eq "t";
	$slurm .= sprintf "module load java/jdk1.8.0_66\n"          if $TNT eq "t";
	$slurm .= sprintf "module load python/3.6-anaconda4.4.0\n"  if $TNT eq "t";
	# parameters
	$slurm .= sprintf "%s/virusgathererHPC-3-assembly.pl %s %s %s -t=%d -a=%s -s=%d -e=%f -i=%d -w=%d -l=%d -k=%d -h=%d -q=%s -p=%d -ol=%d -oh=%d -m=%d -di=%s",
                          $bindir, $family, $project, $assid, $cpusPerTask, $GS_assembler, $GS_seed_size,    $GS_blast_E,     $GS_blast_ide, $GS_blast_ws,
			  $GS_contig_len, $GS_max_itr, $GS_aln_th, $GS_use_qual,  $GS_cap_identity, $GS_cap_overlap, $GS_cap_overhang, $memPerCPUass, $GS_direction;
	$slurm .= " -b"         if ( $useBlastSeed == 1 and $useHitReads  eq "" );
	$slurm .= " -u"         if ( $GS_starter   == 1 and $GS_assembler ne "cap3" );
	$slurm .= " -TNT=$TNT"  if ( $TNT ne "" );
	$slurm .= " -lustre"    if ( $useLustre    == 1 );
	$slurm .= sprintf " -v=%s",  $useHitReads if ( $useHitReads ne "" );
	$slurm .= sprintf " -se=%s", $userSeed    if ( $userSeed    ne "" );
	$slurm .= sprintf " -o=%s",  $userOutdir  if ( $userOutdir  ne "" );
	# create job file
	my $jobfile = "$tmpdir/$assid/$assid-assembly.slurm";
	open (SLURM,">$jobfile");
	print SLURM $slurm."\n";
	close(SLURM);
	# submit command
	my $cmd  = "sbatch";
	   $cmd .= " --dependency=afterok$depjobs" if ( $depjobs ne "" );
	   $cmd .= " $jobfile | awk '{print \$4}'";
	# submit job file and log job ID
	my $jid = `$cmd`;  chomp($jid);
	# delete job file
	#`rm -rf $jobfile`;
	# return job ID
	return( $jid );
}

