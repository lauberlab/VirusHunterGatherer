#!/usr/bin/perl
 
#########################################################
# search for viral sequences in NGS data
# workflow optimzed for usage on a computing server
#
# author: Chris Lauber
#  email: chris.lauber@twincore.de
#########################################################
 
# ------------ #
# load modules
# ------------ #
#use lib "/mnt/twincore/compvironas/software/vhvg/VirusHunterGatherer/lib/perl"; # please make sure to adjust the path
use warnings;
use strict;
# use LWP::Simple qw/get/;
use POSIX qw/floor/;
use Cwd;
use Parallel::ForkManager;

#use lib "/home/lauber/lauber-2015-virusHunter";
#use      virusHGutils;


# --------- #
# parameter
# --------- #
my $hmmerE        = 10;
my $hmmerEsure    = 0.01;
my $hitsTotalHmm  = 2;
#my $refseqDB      = sprintf "%s/db/RefSeq/refseq_protein",        virusHGutils::get_workspace_path( "scratch" );
#my $viralDB       = sprintf "%s/db/RefSeq/viral_genomic",         virusHGutils::get_workspace_path( "scratch" );
#my $refseqNegIDs  = sprintf "%s/db/RefSeq/viral_protein.gi_list", virusHGutils::get_workspace_path( "scratch" );
#my $filterDB      = sprintf "%s/db/RefSeq/filter",                virusHGutils::get_workspace_path( "scratch" );
#my $filterDB      = sprintf "%s/db/Kroepfchen/Kroepfchen_050_entkarpft_reformatted_filtered", virusHGutils::get_workspace_path( "scratch" );
my $id2taxExec    = "1_scripts/queryTax.pl";
my $DBblastEh     = 1e-04;	# cut-off for blast against host sequences
my $DBblastEv     = 1; # 1e-02;	# cut-off for blast against virus sequences
my $cutQth        = 10;	# threshold for cutting low-quality bases in cutadapt
my $assembler     = "cap3";
my $cap3overhang  = 75; # in %
my $cap3overlap   = 20; # in bases
my $cap3ident     = 85; # in %
my $filterE       = 1e-05;
my $filterIdent   = 99;
my $threadsMax    = 12; # don't use more threads than this (to run multiple jobs on very large nodes)
my $transeqbin    = "transeq";
my $transcriptsHg38 = "";
my $genemappingHg38 = "";

my $sraSize;
my $resdir;
my $logdir;
my $proFile;
my $LOGF;

# --------------- #
# usage and input
# --------------- #
if ( $#ARGV != 13 ){
die("
usage: virushunterTWC.pl <parameter>\n
parameters:
\t<family>
\t<fastq_file>
\t<project_id>
\t<run_id>
\t<base_dir>
\t<threads>
\t<filterDB>
\t<refseqDB>
\t<viralDB>
\t<viralRefSeqIDs>
\t<workflow_dir>
\t<filterHg38> [0|1]
\t<debugMode>  [0|1]
\t<enable_filter_3> [0|1]
\n");
}

# -------------------- #
# log program call
# -------------------- #
my $virushuntercall = "virushunterTWC.pl ".join(" ",@ARGV);
print "$virushuntercall\n\n";

# -------------------- #
# get input parameters
# -------------------- #
my $family       = shift;
my $fastqF       = shift;
my $projID       = shift;
my $sraid        = shift;
my $tmpdir       = shift;
my $threads      = shift;
my $filterDB     = shift;
my $refseqDB     = shift;
my $viralDB      = shift;
my $refseqNegIDs = shift;
my $workflowDir  = shift;
my $filterHg38   = shift;
my $debugMode    = shift;
my $enable_filter_3 = shift;


# -------------- #
# initialization
# -------------- #
my ($stime, $rtimeP, $rtimeS, $rtimeF, $rtimeM) = (0,0,0,0,0);
$stime = time();
init();


# ------------ #
# prepare data
# ------------ #
preprocess();


# ----------------------- #
# run HMMer against reads
# ----------------------- #
searchByHMMer();


# ------------------------------------ #
# run Blast against reference proteins
# ------------------------------------ #
searchRefSeq();


# ------------------- #
# compress some files
# ------------------- #
compress();

# --------------------------- #
# total runtime and data size
# --------------------------- #
reportResources();

# ----- #
# clean
# ----- #
cleanup();


# --- #
# END #
# --- #
end();


# --------- #
# functions
# --------- #
# init project (if started the first time)
sub init {
	# create some directories, if not already present
	if( ! -d "$tmpdir/$family" ){
		`mkdir $tmpdir/$family`;
	}
	if( ! -d "$tmpdir/$family/$projID" ){
		`mkdir $tmpdir/$family/$projID`;
	}
	if( ! -d "$tmpdir/$family/$projID/logs" ){
		`mkdir $tmpdir/$family/$projID/logs`;
	}
	if( ! -d "$tmpdir/$family/$projID/logs/virushunter" ){
		`mkdir $tmpdir/$family/$projID/logs/virushunter`;
	}
	if( ! -d "$tmpdir/$family/$projID/results" ){
		`mkdir $tmpdir/$family/$projID/results`;
	}
	if( ! -d "$tmpdir/$family/$projID/results/$sraid" ){
		`mkdir $tmpdir/$family/$projID/results/$sraid`;
	}
	if( ! -d "$tmpdir/$family/$projID/results/$sraid/virushunter" ){
		`mkdir $tmpdir/$family/$projID/results/$sraid/virushunter`;
	}
	# set directories to work with
	$resdir = "$tmpdir/$family/$projID/results/$sraid/virushunter";
	$logdir = "$tmpdir/$family/$projID/logs/virushunter";

	# set databases
	$viralDB   .= "_$family"  if ( -e $viralDB."_$family.nhr"  or  -e $viralDB."_$family.nal" );
	$filterDB  .= "_hg38"     if $filterHg38 != 0;

	# determine size requirements for this run
	$sraSize = `du $fastqF | awk '{print \$1};'`;  chomp( $sraSize );
	$sraSize = $sraSize / 1024 / 1024;

	# set profile file
	$proFile = "$workflowDir/2_profiles/$family-profile.hmm";

	# set paths helper executables
	$id2taxExec = "$workflowDir/$id2taxExec";

	# open log file
	open( $LOGF, ">$logdir/$sraid"."_virushunter.log" ) or die( "Can't write to log file: $!\n" );

	# record name of machine and data directory
	my $hostname = `hostname`; chomp( $hostname );
	printf $LOGF "[virushunter] \t$sraid: running on %s\n", $hostname;
	printf $LOGF "[virushunter] \t$sraid: using data directory %s\n\n", $resdir;

	# limit threads to 32
	$threads = $threadsMax       if ( $threads >= ($threadsMax*2) );
	$threads = int($threads / 2) if ( $threads <= ($threadsMax*2) and $threads > $threadsMax );
	
} # end init

# clean temporary data
sub cleanup {
	# delete temporary files
	if ( $debugMode == 0 ){
		`rm -rf $resdir/$sraid*`;
	}

	# log global statistics
	open(  STATS, ">>$tmpdir/$family/$projID/logs/virushunterTWC-search.txt" ) or die("$0: could not open file: $!\n");
	flock( STATS, 2 ) or die("$0: could not lock file: $!\n");
	printf STATS "%s\t%s\t%.4f\t%d\t%d\n", $sraid, $projID, $sraSize, $rtimeS+$rtimeF, $threads;
	close( STATS ) or die("$0: could not close file: $!\n");

	# log this runs as completed
	open( LOG, ">>$tmpdir/$family/$projID/logs/virushunterTWC-completed.txt" ) or die("$0: could not open file: $!\n");
	flock(LOG, 2) or die("$0: could not lock file: $!\n");
	print LOG $sraid."\n";
	close(LOG) or die("$0: could not close file: $!\n");
	
	# THE END
	printf $LOGF "[virushunter] \t$sraid: done\n\n";
} # end cleanup


# do some last things before ending the scirpt
sub end{
	# put empty contig file if not present (needed for snakemake)
	`touch $resdir/contigs.singlets.fas.gz` if ( ! -e "$resdir/contigs.singlets.fas.gz" );
	# write file showing the search is done
	`echo completed > $resdir/virushunter.done`;
	# close log filehandle
	close( $LOGF );
	# exit with code zero
	exit( 0 );
} # end done


# preprocess data for Blast or HMMer search
sub preprocess {
	printf $LOGF "[virushunter] \t$sraid: transforming to Fasta format\n";
	my $read1id = `zcat $fastqF | head -n 1`;
	my $read2id = `zcat $fastqF | head -n 5 | tail -n 1`;
	chomp($read1id);
	chomp($read2id);
	#$read1id =~ /(.*\.\d+)\.\d .*/;  my $r1id = $1;
	#$read2id =~ /(.*\.\d+)\.\d .*/;  my $r2id = $1;
	my $r1id = "1";
	my $r2id = "2";
	if( $read1id =~ /(.*\.\d+)\.\d .*/ ){
		$r1id = $1;
	}
	if( $read2id =~ /(.*\.\d+)\.\d .*/ ){
		$r2id = $1;
	}
	my $fastp_w   = $threads;
	   $fastp_w   = 16  if $threads > 16;
	my $fastp_cmd = "fastp -i $fastqF -w $fastp_w -j $resdir/$sraid.fastp.json -h $resdir/$sraid.fastp.html";
	if ( $r1id eq $r2id ){
		$fastp_cmd .= " --interleaved_in --stdout > $resdir/$sraid.trim.fastq";
	}else{
		$fastp_cmd .= " -o $resdir/$sraid.trim.fastq";
	}
	`$fastp_cmd 2>/dev/null`;
	`seqtk seq -A $resdir/$sraid.trim.fastq > $resdir/$sraid.fasta 2>/dev/null`;
	`rm $resdir/$sraid.trim.fastq`;
	#exit;
	
	# create 6-frame translation database if needed
	if ( ! -e "$resdir/$sraid-aa-F1.fasta" ){
		printf $LOGF "[virushunter] \t$sraid: creating 6 single-frame translation DBs for HMMer\n";
		my $pm = new Parallel::ForkManager( $threads );
		for ( my $fi=1; $fi<=6; $fi++ ){
			# begin fork
			my $pid   = $pm->start and next;
			# do the work
			my $frame = $fi;
			$frame = ($frame - 3) * -1  if ( $frame > 3 );
			`$transeqbin $resdir/$sraid.fasta $resdir/$sraid-aa-F$fi.fasta -frame=$frame -sformat1=pearson 2>/dev/null`;
			# end fork
			$pm->finish;
		}
		# wait for all forks
		$pm->wait_all_children;
	}
		
	# preprocess runtime for this run
	my $etime = time();
	$rtimeP   = ($etime - $stime);
} # end prepare


# report total runtime and size of SRA dataset
sub reportResources{
	# report runtimes of different parts
	printf $LOGF "\n[virushunter] \t$sraid: runtime preprocess part : %s\n", runtime2hms( $rtimeP );
	printf $LOGF   "[virushunter] \t$sraid: runtime search     part : %s\n", runtime2hms( $rtimeS );
	printf $LOGF   "[virushunter] \t$sraid: runtime filter     part : %s\n", runtime2hms( $rtimeF );
	printf $LOGF   "[virushunter] \t$sraid: total runtime           : %s for %.2f Gb\n", runtime2hms( $rtimeP+$rtimeS+$rtimeF+$rtimeM ), $sraSize;
} # end reportResources

# runtime for h-min-sec format
sub runtime2hms{
	my $rtime = shift;
	my $ho    = floor( $rtime / 3600 );  $rtime -= $ho * 3600;
	my $mi    = floor( $rtime /   60 );  $rtime -= $mi *   60;
	my $hms   = sprintf "%2d h %2d min %2d sec", $ho, $mi, $rtime;
	return( $hms );
} # end runtime2hms


# run and parse HMMer search
sub searchByHMMer {
	# verify that SRA dataset was successfully downloaded
	die( "No data fasta file found!\n" ) if ! -e "$resdir/$sraid.fasta";
	# init new hits and runtime variables
	my %hits  = ();
	# run HMMer
	my $threads_frame = int( $threads / 6 );
	   $threads_frame = 1 if ( $threads_frame < 1 );
	printf $LOGF "[virushunter] \t$sraid: running  6 single-frame HMMer searches using %d threads each\n", $threads_frame;
	my $Ecut = $hmmerE > 10 ? $hmmerE : 10;
	my $targetNum  = `grep -c '>' $resdir/$sraid-aa-F1.fasta`;  chomp( $targetNum );
	   $targetNum *= 6; # set database size to 6 times the number of reads (because we search each reading frame separately) 
	my $pm = new Parallel::ForkManager( 6 );
	for ( my $fi=1; $fi<=6; $fi++ ){
		# begin fgenemappingHg38ork
		my $pid = $pm->start and next;
		# do the work
		my $cmd  = "hmmsearch -E $hmmerE --cpu $threads_frame --tblout $resdir/$sraid-hmmsearch-F$fi.tsv -Z $targetNum $proFile $resdir/$sraid-aa-F$fi.fasta";
		   $cmd .= " > $resdir/$sraid-hmmsearch-F$fi.txt";
		`$cmd`;
		# end fork
		$pm->finish;
	}
	# wait for all forks
	$pm->wait_all_children;
	# combine results of frame-specific hmmsearches
	`cat $resdir/$sraid-hmmsearch-F*.tsv > $resdir/$sraid-hmmsearch.tsv`;
	`cat $resdir/$sraid-hmmsearch-F*.txt > $resdir/$sraid-hmmsearch.txt`;
	# parse HMMer output
	my $evalMin = 1000000;
	open( HMMER, "<$resdir/$sraid-hmmsearch.tsv" );
	while ( my $line = <HMMER> ){
		chomp($line);
		next if $line =~ /^#.*/;
		# get hit read id and E-value
		my @v = split( /\s+/, $line );
		$v[0] =~ /(.*)_\d/;
		my $hid  = $1;
		my $eval = $v[4];
		# save globally best E-value
		if ( $eval < $evalMin ){  $evalMin = $eval; }
		# save hit if E-value below threshold
		if ( $eval < $hmmerE ){
			$hits{$hid}{'score'}    = $v[5];
			$hits{$hid}{'eval'}     = $eval;
			$hits{$hid}{'qProfile'} = $v[2];
		}
	}
	close(HMMER);
	# determine if this SRA is a hits and save HMMer output for it
	my $numhits = scalar keys %hits;
	my $iamhit = "false";
	if ( ( $numhits >= $hitsTotalHmm ) or
	     ( $evalMin < $hmmerEsure ) ){
		$iamhit = "true";
		# save full HMMer output
		`mv $resdir/$sraid-hmmsearch.txt $resdir/hmmsearch.out`;
		# save hits separately
		open( HITS, ">$resdir/hmmsearch-hits.tsv" );
		foreach ( sort {$hits{$a}{'eval'} <=> $hits{$b}{'eval'}} keys %hits ){
			printf HITS "%s\t%s\t%s\t%s\n", $_, $hits{$_}{'qProfile'},
			                                $hits{$_}{'score'}, $hits{$_}{'eval'};
		}
		close(HITS);
	}else{
		`rm -rf $resdir/$sraid.fasta`;
		`rm -rf $resdir/*` if ( $debugMode == 0 );
	}
	
	# search runtime for this run
	my $etime = time();
	$rtimeS   = ($etime - $stime - $rtimeP);
	
	# report final info for this run
	printf $LOGF "[virushunter] \t$sraid: hmmsearch hits:            %d\n",   $numhits;
	printf $LOGF "[virushunter] \t$sraid: best E-value:              %.1e\n", $evalMin;
	printf $LOGF "[virushunter] \t$sraid: detected as hit:           %s\n",   $iamhit;
} # end searchByHMMer


# step-wise filtering of hit reads using RefSeq
sub searchRefSeq{
	# only run for SRA experiments that are hits
	return 0 if ! -d "$resdir";
	# verify that SRA dataset was successfully prepared
	return 0 if ! -e "$resdir/$sraid.fasta";
	printf $LOGF "\n[virushunter] \t$sraid: comparing hits against host and viral reference sequences\n";
        # read hits and info from search result file
        my $rawHitsFile = "$resdir/hmmsearch-hits.tsv";
	my %readHits = ();
	open( RH, "<$rawHitsFile" ) or die( "Can't open file '$rawHitsFile': $!\n" );
	while ( my $line = <RH> ){
		chomp($line);  next if $line eq "";
		my @v = split( /\t/, $line );
		$readHits{ $v[0] } = ();
		$readHits{ $v[0] }{ 'query' } = $v[1];
		$readHits{ $v[0] }{ 'E' }     = $v[3];
	}
	close(RH);
	# extract reads
	my $rawHitsIDs = "$resdir/$sraid-HitReadsAll.ids";
	my $rawHitsFas = "$resdir/$sraid-HitReadsAll.fasta";
	open( RHIDS, ">$rawHitsIDs" ) or die( "Can't open file '$rawHitsIDs': $!\n" );
	foreach ( sort keys %readHits ){  printf RHIDS "%s\n", $_; }
	close( RHIDS);
	my $cmd1 = "grep -F --no-group-separator -A 1 -f $rawHitsIDs $resdir/$sraid.fasta > $rawHitsFas";
	`$cmd1`;
	# dereplicate reads
	printf $LOGF "[virushunter] \t$sraid: dereplicating initial hits with vsearch\n";
	my $rawHitsFasDerep = "$resdir/$sraid-HitReadsAll-dereplicated.fasta";
	my $cmdDerep = "vsearch --derep_fulllength $rawHitsFas --notrunclabels --output $rawHitsFasDerep --quiet";
	`$cmdDerep`;
	my $nHitsLeft = `grep -c ">" $rawHitsFasDerep`;
    chomp($nHitsLeft);
	printf $LOGF "[virushunter] \t$sraid: $nHitsLeft initial hits left\n";
	# assemble reads to contigs
	printf $LOGF "[virushunter] \t$sraid: assembling initial hits with $assembler\n";
	my $assLog      = $rawHitsFasDerep.".cap.log";
	my $assContigs  = $rawHitsFasDerep.".cap.contigs";
	my $assSinglets = $rawHitsFasDerep.".cap.singlets";
	my $cmd2 = "$assembler $rawHitsFasDerep -h $cap3overhang -o $cap3overlap > $assLog";		
	`$cmd2`;
	# read contigs and reformat
	my %contigs  = ();
	my $contigID = "";
	open( CNTG, "<$assContigs" ) or die ( "Can't open file '$assContigs': $!\n" );
	while ( my $line = <CNTG> ){
		chomp($line);  next if $line eq "";
		if ( $line =~ />(.*)/ ){
			$contigID = $sraid.".".$1.".1";
			$contigs{ $contigID } = {};
			$contigs{ $contigID }->{ 'seq' }       = "";
			$contigs{ $contigID }->{ 'reads' }     = [];
			$contigs{ $contigID }->{ 'E_refseq' }  = 1e6;		# init here, to be able to switch filters   
			$contigs{ $contigID }->{ 'inner_contigID' } = $contigID;	# init here, to be able to switch filters   
		}else{
			$contigs{ $contigID }->{ 'seq' }      .= $line;
		}
	}
	close(CNTG);
	# get read-contig mapping
	open( CAP, "<$assLog" ) or die( "Can't open file '$assLog': $!\n" );
	$contigID = "";
	while ( my $line = <CAP> ){
		chomp($line);  next if $line eq "";
		last if $line eq "DETAILED DISPLAY OF CONTIGS";
		if ( $line =~ /\*+ Contig (\d+) \*+/ ){
			$contigID = $sraid.".Contig".$1.".1";
		}elsif ( $contigID ne "" ){
			$line =~ s/^\s+//; $line =~ s/ .*//;  $line =~ s/\-$//;  $line =~ s/\+$//;
			push( @{ $contigs{ $contigID }->{ 'reads' } }, $line );
		}		
	}
	close(CAP);
	# read singlets and combine with contigs
	open( SNGLT, "<$assSinglets" ) or die ( "Can't open file '$assSinglets': $!\n" );
	while ( my $line = <SNGLT> ){
		chomp($line);  next if $line eq "";
		if ( $line =~ />(.*)/ ){
			$contigID = $1;  $contigID =~ s/ .*//;
			$contigs{ $contigID }->{ 'seq' }       = "";
			$contigs{ $contigID }->{ 'reads' }     = [];
			$contigs{ $contigID }->{ 'E_refseq' }  = 1e6;		# init here, to be able to switch filters   
			$contigs{ $contigID }->{ 'inner_contigID' } = $contigID;	# init here, to be able to switch filters   
			push( @{ $contigs{ $contigID }->{ 'reads' } }, $contigID );
		}else{
			$contigs{ $contigID }->{ 'seq' }      .= $line;
		}
	}
	close(SNGLT);
	# write reformatted contigs and singlets to file
	my $contigIDsRaw = "$resdir/contigs.singlets-unfiltered.txt";
	open( CONTSING, ">$contigIDsRaw" ) or die ( "Can't open file '$contigIDsRaw': $!\n" );
	foreach ( sort keys %contigs ){
		printf CONTSING "%s\n", $_;
	}
	close(CONTSING);
	# run filter 1 - blastx against custom contaminant DB
	my $contigIDsFilter1  = "$resdir/contigs.singlets-filter1.txt";
	my $contigSeqsFilter1 = "$resdir/$sraid-HitReadsAll-trimmed.fasta.cap.contigs.singlets-filter1.fas";
	filter1( $contigIDsRaw,     $contigIDsFilter1, $contigSeqsFilter1, \%contigs );

	my $contigIDsFilter3  = "$resdir/contigs.singlets-filter3.txt";
	my $contigSeqsFilter3 = "$resdir/$sraid-HitReadsAll-trimmed.fasta.cap.contigs.singlets-filter3.fas";


	my $contigIDsFilter2  = "$resdir/contigs.singlets-filter2.txt";
	my $contigSeqsFilter2 = "$resdir/$sraid-HitReadsAll-trimmed.fasta.cap.contigs.singlets-filter2.fas";

	if ($enable_filter_3) {
    		# If filter 3 is enabled, run filter 3 next
		# run filter 3 - blastx against non-viral refseq_protein
		# sergej: added to skip filter
		printf $LOGF "[virushunter]   Filter 3 is enabled.\n";
	
    		filter3( $contigIDsFilter1, $contigIDsFilter3, $contigSeqsFilter3, \%contigs );

    		# Then run filter 2 using the output of filter 3
		# run filter 2 - tblastx against viral-genomic
    		filter2( $contigIDsFilter3, $contigIDsFilter2, $contigSeqsFilter2, \%contigs );
	} else {
    		# If filter 3 is not enabled, skip it and run filter 2 using the output of filter 1
		# run filter 2 - tblastx against viral-genomic
		printf $LOGF "[virushunter]   Filter 3 is not enabled.\n";
    		filter2( $contigIDsFilter1, $contigIDsFilter2, $contigSeqsFilter2, \%contigs );
	}
	# extract remaining hits to file
	my $contigSeqsFinal   = "$resdir/contigs.singlets.fas";
	my $contigReadsFinal  = "$resdir/contigs.singlets.reads.tsv";
	my @idsFinal = ();
	open( SOUT, ">$contigSeqsFinal" )  or die ( "Can't open file '$contigSeqsFinal': $!\n" );	
	open( ROUT, ">$contigReadsFinal" ) or die ( "Can't open file '$contigReadsFinal': $!\n" );	
	open(  IIN, "<$contigIDsFilter2" ) or die ( "Can't open file '$contigIDsFilter2': $!\n" );
	while ( my $line = <IIN> ){
		chomp($line);  next if $line eq "";
		push( @idsFinal, $line );
		printf SOUT ">%s\n%s\n", $line, $contigs{ $line }->{ 'seq' };
		foreach ( @{ $contigs{ $line }->{ 'reads' } } ){
			printf ROUT "%s\t%s\n", $line, $_;
		}
	}	
	close( IIN );
	close(ROUT );
	close(SOUT );
	# remove this SRA from hit directory if no hits remain after filtering
	if ( scalar @idsFinal == 0 ){
		# filter runtime for this run
		my $etime = time();
		$rtimeF   = ($etime - $stime - $rtimeP - $rtimeS);
		# end here
		`rm -rf $resdir/*` if ( $debugMode == 0 );
		return 0;
	}
	# summarize final hits
	my $res = {};
	my $contig2query = {};
	foreach my $sid ( @idsFinal ){
		# sort by viral refseq subject that was hit
		my $vid = $contigs{ $sid }->{ 'id_refseq' } // ""; # sergej: added // incase where no info can be fetched
		if ( ! exists $res->{ $vid } ){
			$res->{ $vid } = {};
			$res->{ $vid }->{ 'init_best_E' }       = 1e6;
			$res->{ $vid }->{ 'init_best_query' }   = "";
			$res->{ $vid }->{ 'init_num_hits' }     = 0;
			$res->{ $vid }->{ 'refseq_contigs' }    = 0;
			$res->{ $vid }->{ 'refseq_best_E' }     = 1e6;
			$res->{ $vid }->{ 'refseq_best_Ident' } = 0;
			# it will look like 'NC_006430.1|Mammalian orthorubulavirus 5'
			$res->{ $vid }->{ 'refseq_subject' }    = sprintf "%s|%s", $vid, $contigs{ $sid }->{ 'title_refseq' }// ""; # sergej: added // for cases where information cannot be fetched to prevent warnings.
		}
		# best E, query, and number of hit reads of initial search
		$contig2query->{ $sid } = {};
		$contig2query->{ $sid }->{ 'E' }     = 1e6;
		$contig2query->{ $sid }->{ 'query' } = "";
		foreach my $rid ( @{ $contigs{ $sid }->{ 'reads' } } ){
			if ( $readHits{ $rid }{ 'E' } < $res->{ $vid }->{ 'init_best_E' } ){
				$res->{ $vid }->{ 'init_best_E' }     = $readHits{ $rid }{ 'E' };
				$res->{ $vid }->{ 'init_best_query' } = $readHits{ $rid }{ 'query' };
			}
			$res->{ $vid }->{ 'init_num_hits' }++;
			if ( $readHits{ $rid }{ 'E' } < $contig2query->{ $sid }->{ 'E' } ){
				$contig2query->{ $sid }->{ 'E' }     = $readHits{ $rid }{ 'E' };
				$contig2query->{ $sid }->{ 'query' } = $readHits{ $rid }{ 'query' };
			}
		}
		# best E, subject, number of reads and contigs, and subject name of RefSeq search
		if ( $contigs{ $sid }->{ 'E_refseq' } < $res->{ $vid }->{ 'refseq_best_E' } ){
			$res->{ $vid }->{ 'refseq_best_E' }     = $contigs{ $sid }->{ 'E_refseq' };
			$res->{ $vid }->{ 'refseq_best_Ident' } = $contigs{ $sid }->{ 'ident_refseq' };
			$res->{ $vid }->{ 'refseq_best_lens' }  = $contigs{ $sid }->{ 'lens_refseq' };
		}
		$res->{ $vid }->{ 'refseq_contigs' }++;
	}
	# annotate final hits with taxonomy info
	printf $LOGF "[virushunter] \t$sraid: adding taxonomy annotation for final viral reference sequence hits\n";	
	my $idFile = "$resdir/$sraid-FinalHits-ids.txt";
	open( ID, ">$idFile" ) or die( "Can't open file '$idFile': $!\n" );
	foreach my $vid ( keys %{$res} ) {
		my $refseq_subject = $res->{ $vid }->{ 'refseq_subject' };
		# refseq_subject looks like 'NC_006430.1|Mammalian orthorubulavirus 5'
		my ($id_refseq) = $refseq_subject =~ /^([^|]+)/;
		$id_refseq //= ''; # sergej: added to supress warning for when undefined
		print ID "$id_refseq\n";
	}
	close(ID);
	my $cmdTax = "$id2taxExec $idFile";
	open( TAX, "$cmdTax |" );
	while ( my $line = <TAX> ){
		chomp($line);  next if $line eq "";
		my @v = split( /\t/, $line );
		$res->{ $v[0] }->{ 'refseq_taxonomy' } = $v[1];
	}
	close(TAX);
	# write contig to initial_query mapping to file
	my $contig2queryFile = "$resdir/contigs.singlets.query.tsv";
	open( C2Q, ">$contig2queryFile" ) or die( "Can't open file '$contig2queryFile': $!\n" );
	foreach my $sid ( sort { $contig2query->{ $a }->{ 'query' } cmp $contig2query->{ $b }->{ 'query' } } keys %{ $contig2query } ){
		printf C2Q "%s\t%s\t%.2e\n", $sid, $contig2query->{ $sid }->{ 'query' }, $contig2query->{ $sid }->{ 'E' } ;
	}
	close(C2Q);	
	# write result summary to file
	my $finalResFile = "$resdir/final.hits.tsv";
	open( RES, ">$finalResFile" ) or die( "Can't open file '$finalResFile': $!\n" );
	foreach my $vid ( sort { $res->{$a}->{'init_best_E'} <=> $res->{$b}->{'init_best_E'} } keys %{$res} ){
		printf RES "%d\t%.2e\t%s\t%.2e\t%.1f\t%s\t%d\t%s\t", $res->{$vid}->{'init_num_hits'},  $res->{$vid}->{'init_best_E'},       $res->{$vid}->{'init_best_query'},
							  	     $res->{$vid}->{'refseq_best_E'},  $res->{$vid}->{'refseq_best_Ident'}, $res->{$vid}->{'refseq_best_lens'},
								     $res->{$vid}->{'refseq_contigs'}, $res->{$vid}->{'refseq_subject'};
		printf RES "%s",                                     $res->{$vid}->{'refseq_taxonomy'}  if exists $res->{$vid}->{'refseq_taxonomy'};
		printf RES "\n";
	}
	close(RES);
	# filter runtime for this run
	my $etime = time();
	$rtimeF   = ($etime - $stime - $rtimeP - $rtimeS);
}

# filter - blast against custom contaminant DB
sub filter1{
	my $idsin  = shift;
	my $idsout = shift;
	my $fasout = shift;
	my $cntgs  = shift;
	#printf $LOGF "[virushunter] \t$sraid: filtering - tblastx against custom contaminant DB\n";
	printf $LOGF "[virushunter] \t$sraid: filtering -  blastn against custom contaminant DB\n";
	# read query sequence IDs from file
	open( IDSIN, "<$idsin" ) or die ( "Can't open file '$idsin': $!\n" );
	my %ids = ();
	while ( my $line = <IDSIN> ){
		chomp($line);  next if $line eq "";  $ids{ $line } = 1;
	}
	close(IDSIN);
	# write query sequences to file
	open( FASOUT, ">$fasout" ) or die ( "Can't open file '$fasout': $!\n" );
	foreach ( sort keys %ids ){
		printf FASOUT ">%s\n%s\n", $_, $$cntgs{ $_ }->{ 'seq' };
	}
	close(FASOUT);
	# blast
	#my $blastfile = $fasout."-tblastx.tsv";
	my $blastfile = $fasout."-blastn.tsv";
	#my $cmd = "tblastx -db $16543.hunter_testfilterDB -query $fasout -outfmt '6 qseqid sseqid pident evalue' -evalue $DBblastEh -num_threads $threads > $blastfile";
	my $cmd = "blastn -db $filterDB -query $fasout -outfmt '6 qseqid sseqid pident evalue' -evalue $DBblastEh -num_threads $threads > $blastfile";
	`$cmd 2>/dev/null`;
	# get hits
	open( BLAST, "<$blastfile" ) or die( "Can't open file '$blastfile': $!\n" );
	while ( my $line = <BLAST> ){
		chomp($line);  next if $line eq "";
		my @v = split( /\t/, $line );
		delete $ids{ $v[0] };
	}
	close(BLAST);
	# write positive hit sequences to file
	open( IDSOUT, ">$idsout" ) or die( "Can't open file '$idsout': $!\n" );
	foreach ( sort keys %ids ){
		printf IDSOUT "%s\n", $_;
	}
	close(IDSOUT);
	# report remaining number of sequences
	printf $LOGF "[virushunter] \t$sraid: filtering -  %d hits remaining\n", scalar keys %ids;
} # end filter1

# filter - tblastx against viral_genomic
sub filter2{
	my $idsin  = shift;
	my $idsout = shift;
	my $fasout = shift;
	my $cntgs  = shift;
	printf $LOGF "[virushunter] \t$sraid: filtering - tblastx against viral_genomic\n";
	# read query sequence IDs from file
	open( IDSIN, "<$idsin" ) or die ( "Can't open file '$idsin': $!\n" );
	my %ids = ();
	while ( my $line = <IDSIN> ){
		chomp($line);  next if $line eq "";  $ids{ $line } = 1;
	}
	close(IDSIN);
	# id2taxExec write query sequences to file
	open( FASOUT, ">$fasout" ) or die ( "Can't open file '$fasout': $!\n" );
	foreach ( sort keys %ids ){
		printf FASOUT ">%s\n%s\n", $_, $$cntgs{ $_ }->{ 'seq' };
	}
	close(FASOUT);
	# blast
	my $blastfile = $fasout."-tblastx.tsv";
	my $cmd = "tblastx -db $viralDB -query $fasout -outfmt '6 qseqid stitle pident evalue length qlen sseqid' -max_hsps 1 -max_target_seqs 1 -evalue $DBblastEv -num_threads $threads > $blastfile";
	`$cmd 2>/dev/null`;
	# get hits
	my %idsO = ();
	open( BLAST, "<$blastfile" ) or die( "Can't open file '$blastfile': $!\n" );
	while ( my $line = <BLAST> ){
		chomp($line);  next if $line eq "";
		my @v = split( /\t/, $line );
		$idsO{   $v[0] } = 1;
		$$cntgs{ $v[0] }->{ 'E_refseq' }      = $v[3];
		$$cntgs{ $v[0] }->{ 'ident_refseq' }  = $v[2];
		# sseqid looks like 'ref|NC_006430.1|' - extact the actual id between pipes
		( $$cntgs{ $v[0] }->{ 'id_refseq' } ) = $v[6] =~ /\|(.+)\|\z/;
        # if empty, let's save the whole $v[6]
		# Why? It could be just sequence id, not surrounded by pipes (this is the case, when
		# additional genomes are included manually into the RefSeq viral_genomes database).
		$$cntgs{ $v[0] }->{ 'id_refseq' } ||= $v[6];
		$$cntgs{ $v[0] }->{ 'title_refseq' }  = $v[1];
		$$cntgs{ $v[0] }->{ 'lens_refseq' }   = sprintf "%d / %d", $v[4]*3, $v[5];
	}
	close(BLAST);
	# save hits
	`mv $blastfile $resdir/tblastx-viral_genomic-hits.tsv`;
	# write positive hit sequences to file
	open( IDSOUT, ">$idsout" ) or die( "Can't open file '$idsout': $!\n" );
	foreach ( sort keys %idsO ){
		printf IDSOUT "%s\n", $_;
	}
	close(IDSOUT);
	# report remaining number of sequences
	printf $LOGF "[virushunter] \t$sraid: filtering -  %d hits remaining\n", scalar keys %idsO;
} # end filter2

# filter - blastx against non-viral refseq_protein
sub filter3{
	my $idsin  = shift;
	my $idsout = shift;
	my $fasout = shift;
	my $cntgs  = shift;
	printf $LOGF "[virushunter] \t$sraid: filtering -  blastx against non-viral refseq_protein\n";
	# read query sequence IDs from file
	open( IDSIN, "<$idsin" ) or die ( "Can't open file '$idsin': $!\n" );
	my %ids = ();
	while ( my $line = <IDSIN> ){
		chomp($line);  next if $line eq "";  $ids{ $line } = 1;
	}
	close(IDSIN);
	# only analyze contig with best E value per refseq subject
	my %ids2 = ();
	foreach ( keys %{$cntgs} ){
		my $vid = $$cntgs{ $_ }->{ 'inner_contigID' };
		if ( ! exists $ids2{ $vid } ){
			$ids2{ $vid } = ();
			$ids2{ $vid }{ 'contigID' } = $_;
			$ids2{ $vid }{ 'Evalue' }   = $$cntgs{ $_ }->{ 'E_refseq' };
		}elsif( $$cntgs{ $_ }->{ 'E_refseq' } < $ids2{ $vid }{ 'Evalue' } ){
			$ids2{ $vid }{ 'contigID' } = $_;
			$ids2{ $vid }{ 'Evalue' }   = $$cntgs{ $_ }->{ 'E_refseq' };
		}
	}
	# write query sequences to file
	open( FASOUT, ">$fasout" ) or die ( "Can't open file '$fasout': $!\n" );
	foreach ( sort keys %ids2 ){
		my $cid = $ids2{ $_ }{ 'contigID' };
		printf FASOUT ">%s\n%s\n", $cid, $$cntgs{ $cid }->{ 'seq' };
	}
	close(FASOUT);
	# blast
	my $blastfile = $fasout."-blastx.tsv";
	my $cmd  = "blastx -db $refseqDB -query $fasout -outfmt '6 qseqid sacc stitle pident evalue' -max_hsps 1 -max_target_seqs 1 -evalue $DBblastEh -num_threads $threads ";
	  #$cmd .= "-negative_gilist $refseqNegGIs > $blastfile";
	   $cmd .= "-negative_seqidlist $refseqNegIDs > $blastfile.0";
	`$cmd 2>/dev/null`;
	`grep -v 'hypothetical protein' $blastfile.0 > $blastfile`;
	# get hits
	open( BLAST, "<$blastfile" ) or die( "Can't open file '$blastfile': $!\n" );
	while ( my $line = <BLAST> ){
		chomp($line);  next if $line eq "";
		my @v = split( /\t/, $line );
		next if $v[2] =~ /.*virus.*/i;
		my $cid = $v[0];
		my $vid = $$cntgs{ $cid }->{ 'inner_contigID' };
		# delete all contigs that hit the respective refseq subject
		foreach ( keys %{$cntgs} ){
			if ( $$cntgs{ $_ }->{ 'inner_contigID' } eq $vid ){
				delete $ids{ $_ };
			}
		}
	}
	close(BLAST);
	# write positive hit sequences to file
	open( IDSOUT, ">$idsout" ) or die( "Can't open file '$idsout': $!\n" );
	foreach ( sort keys %ids ){
		printf IDSOUT "%s\n", $_;
	}
	close(IDSOUT);
	# report remaining number of sequences
	printf $LOGF "[virushunter] \t$sraid: filtering -  %d hits remaining\n", scalar keys %ids;
} # end filter3

# trim reads using autoadapt
# ! currently not used n this workflow !
sub trimReads{
	my $infile  = shift;
	my $outfile = shift;
	printf $LOGF "[virushunter] \t$sraid: trimming adapter sequences and low-quality bases for initial hits\n";
	# prepare formatted read IDs
	my $readIDfile = "$infile.fqids";
	`grep '>' $infile | sed 's/>/@/' > $readIDfile`;
	# temporary extract fastq reads
	my $fastqFile1 = "$infile";
	   $fastqFile1 =~ s/\.fasta/\.fastq/;
	`fastq-dump -B --split-spot --skip-technical --readids -O $tmpdir $tmpdir/$sraid.sra`;
	`grep -F -f $readIDfile --no-group-separator -A 3 $tmpdir/$sraid.fastq > $fastqFile1`;
	`rm $tmpdir/$sraid.fastq`;
	# run autoadapt
	my $fastqFile2   = "$fastqFile1";
	   $fastqFile2   =~ s/\.fastq/-trimmed\.fastq/;
	my $autoadaptLog = "$infile.autoadapt.log";
	my $wdir = cwd;
	chdir( $tmpdir );
	`autoadapt --threads=$threads --quality-cutoff=$cutQth $fastqFile1 $fastqFile2 > $autoadaptLog`;
	chdir( $wdir );
	# convert fastq to fasta
	`seqtk seq -a $fastqFile2 > $outfile`;
} # end trimReads


# map all reads against reference transcriptome
# ! not yet included into this workflow !
sub mapRefGenes{
	# only run for SRA experiments that are hits
	return 0 if ! -d "$resdir";
	# verify that SRA dataset was successfully downloaded
	return 0 if ! -e "$resdir/$sraid.fasta";
	printf $LOGF "\n[virushunter] \t$sraid: mapping reads against reference transcripts\n";
	# quantify via salmon
	my $cmd  = "salmon quant -l A -i $transcriptsHg38 -r $tmpdir/$sraid.fasta -p $threads";
	   $cmd .= " -g $genemappingHg38 -o $resdir/$sraid-hg38-salmon > $resdir/$sraid-hg38-salmon.log 2>&1";
	`$cmd`;
	# save result file
	#`cp $tmpdir/$sraid-hg38-salmon/quant.genes.sf $tmpdir/$projID/$sraid/hg38.counts`;
	`awk '{print \$1,\$4,\$5}' $resdir/$sraid-hg38-salmon/quant.genes.sf > $resdir/hg38.counts`;
	# mapping runtime for this run
	my $etime = time();
	$rtimeM   = ($etime - $stime - $rtimeP - $rtimeS - $rtimeF);
} # end mapRefGenes


# function to compress some result files to save disk space
sub compress{
	# Blast/HMMer output
	if ( -e "$resdir/hmmsearch.out" ){                  `gzip $resdir/hmmsearch.out; gzip $resdir/hmmsearch-hits.tsv`; }
	if ( -e "$resdir/tblastx-viral_genomic-hits.tsv" ){ `gzip $resdir/tblastx-viral_genomic-hits.tsv`; }
	# other result files
	if ( -e "$resdir/contigs.singlets.fas" ){        `gzip $resdir/contigs.singlets.fas`; }
	if ( -e "$resdir/contigs.singlets.reads.tsv" ){  `gzip $resdir/contigs.singlets.reads.tsv`; }
	if ( -e "$resdir/contigs.singlets.query.tsv" ){  `gzip $resdir/contigs.singlets.query.tsv`; }
	if ( -e "$resdir/hg38.counts" ){                 `gzip $resdir/hg38.counts`; }
} # end compress

