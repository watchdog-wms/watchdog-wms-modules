#!/usr/bin/env perl

use strict;

use File::Basename qw(dirname);
use File::Path qw(mkpath);
use lib dirname (__FILE__).'/../../core_lib/perl';

use IO::Handle;
use Watchdog::ExitCode qw(exitCode);
use Getopt::Long::Descriptive;

# the script adds sequences to SAM files.
# depending on the size of the fastq file the script needs a lot of RAM because all sequence identifiers are load into a HashMap!

my ($opt, $usage) = describe_options(
'addSequences2Sam.pl %o',
[ 'sam|s=s',		"path to the SAM file", { required => 1 }],
[ 'fastq|f=s',		"path to the FASTQ file, multiple files separated with ','", { required => 1 }],
[ 'output|o=s',		"path to the output SAM file in which the sequences are added", { required => 1 }],
[ 'unmapped|u=s',	"[optional] path to a FASTQ file in which the unmapped sequences will be written to; exclusive with --preread flag"],
[ 'noquality',		"[optional] does not add the read quality values; default: disabled"],
[ 'update', 		"[optional] overrides already existing output files; default: disabled"],
[ 'preread', 		"[optional] does only index reads stored in the FASTQ file that are part of the SAM file; exclusive with --unmapped parameter; default: disabled"],
[ 'returnfilepath|r=s', "path to the return variables file"], 
[ 'help',		"print usage message and exit" ]
);
print($usage->text), exitCode("EXIT_OK") if $opt->help;

local $| = 1; # flush STDOUT
(my $filenameSAM = $opt->sam) =~ s/['"]//g;
(my $filenameFastq = $opt->fastq) =~ s/['"]//g;
(my $filenameOUT = $opt->output) =~ s/['"]//g;
(my $filenameUnmapped = $opt->unmapped) =~ s/['"]//g;
(my $quality = !$opt->noquality) =~ s/['"]//g; 
(my $update = $opt->update) =~ s/['"]//g;
(my $preread = $opt->preread) =~ s/['"]//g;
(my $returnFilePath = $opt->returnfilepath) =~ s/['"]//g;

if(!defined $filenameUnmapped || $filenameUnmapped eq "") {
	$filenameUnmapped = "-1";
}

# test, if exclusive flags are used
if($filenameUnmapped ne "-1" && $preread) {
	print "Flag --preread and parameter --unmapped are exclusive.\n";
	exitCode("EXIT_INVALID_ARGUMENTS");
}

open(my $fh_sam, '<:encoding(UTF-8)', $filenameSAM) or do { print "Could not open sam file '$filenameSAM'!\n"; exitCode("EXIT_MISSING_INPUT_FILES"); };
my $fh_out;
my @fastqFiles=split(",", $filenameFastq);

if(!-d $filenameOUT) {
	mkpath(dirname($filenameOUT)); 
}

# test, if file is already there 
if(-e $filenameOUT) {
	# delete the file
	if($update) {
		open($fh_out, '>:encoding(UTF-8)', $filenameOUT) or do { print "Could not open output file '$filenameOUT'!\n"; exitCode("EXIT_WRITING_FAILED"); };
		close($fh_out);
		unlink $filenameOUT;
	}
	else {
		print "Old file there, skipped processing. Call with --update flag in order to overwrite the old file.\n";
		if($filenameUnmapped eq "-1") {
			$filenameUnmapped = "";
		}
		writeReturnValues($filenameOUT, $filenameUnmapped);
		exitCode("EXIT_OK");
	}
}

my $fh_unmapped;
if($filenameUnmapped ne "-1") {
	if(!-d $filenameUnmapped) {
		mkpath(dirname($filenameUnmapped));
	}
	open($fh_unmapped, '>:encoding(UTF-8)', $filenameUnmapped) or do { print "Could not open output file '$filenameUnmapped'!\n"; exitCode("EXIT_WRITING_FAILED"); };
	# is ok, delete the file and reopen it 
	close($fh_unmapped);
	unlink $filenameUnmapped;
}
print "SAM file: $filenameSAM\n";
print "Output file: $filenameOUT\n";
if($filenameUnmapped ne "-1") {
	print "Unmapped FASTQ file: $filenameUnmapped\n";
}

# index IDs based on SAM file
my %index;
if($preread) {
	print "started indexing of mapped reads part of the SAM file...\n";
	my $id;
	my $flag;
	my @tmp;
	my $row;
	while ($row = <$fh_sam>) {
		chomp $row;
		if($row !~ /^@/)  {  
			@tmp = split("\t", $row); 
			if($#tmp >= 10) {
				$id = $tmp[0];
				$flag = $tmp[1];
				$id = getReadIDWithSuffix($id, $flag);

				# save that id in the index hash to indicate that it should be indexed afterwards
				$index{$id} = 1;
			}
		}
	}
	# reset file pointer
	seek($fh_sam, 0, 0);
	print scalar(keys(%index))." reads are part of the SAM file and will be indexed now...\n";
}

my @mappingStatus = ();
my @fileHandles = ();
my $fcount = 0;
my $i = 0;
print "started indexing of reads...\n";
my $fC = 0;
foreach my $fFile (@fastqFiles) {
	$fC++;
	open(my $fh_fastq, '<:encoding(UTF-8)', $fFile) or do { print "Could not open fastq file '$fFile'!\n"; exitCode("EXIT_MISSING_INPUT_FILES"); };
	$fileHandles[$fcount] = $fh_fastq;

	print "FastQ file: $fFile\n";
	my $pos = 0;
	my $id;
	my $row;
	my $l; 
	my $lSize = 0;
	my $lc = 0;

	# check how big a file ending is (\n or \r\n)
	if(defined($row = <$fh_fastq>)) {
		chomp $row;
		$l = length($row);
		$lSize = tell($fh_fastq) - $l;
	}

	# create index for random file access
	my $indexNumber;
	seek($fh_fastq, 0, 0);
	while ($row = <$fh_fastq>) {
		$lc++;
		# get only ID lines
		next unless ($lc % 4 == 1); 
		chomp $row;

		# check, if fastq ID
		if($row =~ /^@/)  {

			# save position of sequence for id
			$id = substr($row, 1);

			# replace it with the old format
			if($id =~ /(.+(\/[12])?)\s+([0-9]):[NY]:[0-9]+:#?([ATCGN]{6})?/) {
				$id = $1;
				$indexNumber = $3;
				if($indexNumber > 2) { # fix for RNA-seq Raji Eick
					$indexNumber = 2;
				}
				if($1 !~ /\/[12]$/) {
					$id = "$id/$indexNumber"
				}
			}
			# assume $fC as index
			if($id !~ m!/[12]$!) {
				$id = "$id/$fC"
			}
		
			# save position if not in preread mode or if read with that ID should be saved
			if(!$preread || exists $index{$id}) {
				$pos = tell($fh_fastq);
				$mappingStatus[$i] = [$pos-(length($id)+1)-$lSize, 0, $fcount];
				$index{$id} = [$pos, $i, $fcount];
			}
			$i++;

			# status update
			if($i % 250000 == 0) { 
				print "processed ".($i)." reads (last ID: $id)...\n";
			}
		}
		else {
			print "No FastQ ID was found in line '$lc': '$row'!\n";
			exitCode("EXIT_MISFORMATED_INPUT");
		}
	}
	print "finished indexing of ".($lc/4)." reads for file '$fFile'!\n";
	$fcount++;
}

open($fh_out, '>:encoding(UTF-8)', $filenameOUT) or do { print "Could not open output file '$filenameOUT'!\n"; exitCode("EXIT_WRITING_FAILED"); };
print "started parsing of SAM file...\n";
# read sam file
my @tmp;
my @tmp2;
my $seq;
my $qual;
my $cigar;
my $notFound = 0;
my $usedID = 0;
my $flag;
my $cstart = 0;
my $cend = 0;
my $reverse = 0;
my $i = 0;
my $row;
my $id;
my $fh_fastq;
my @data;
my $cc = 0;
while ($row = <$fh_sam>) {
$cc++;
	chomp $row;
	# output header unprocessed
	if($row =~ /^@/)  {
		print $fh_out $row."\n";
	}
	else {  
		@tmp = split("\t", $row); 
		if($#tmp >= 10) {
			$id = $tmp[0];
			$flag = $tmp[1];
			$id = getReadIDWithSuffix($id, $flag);

			# get sequence out of fastq file
			if(exists $index{$id}) {
				@data = @{$index{$id}};
				# select correct file handle
				$fh_fastq = $fileHandles[$data[2]];
				# jump to correct position
				seek($fh_fastq, $data[0], 0);
				$seq = <$fh_fastq>;
				chomp $seq;

				# reverse complement if the flag is set
				if($flag & 16) {
					$seq = revCom($seq);
					$reverse = 1;
				}
				
				# check for hard clipping
				$cigar = $tmp[5];
				if($cigar =~ /^([0-9]+)H/) { #hard clipping at start
					$seq = substr($seq, $1);
					$cstart = $1;
				}
				if($cigar =~ /([0-9]+)H$/) { #hard clipping at end
					$seq = substr($seq, 0, -$1);
					$cend = -$1;
				}

				$tmp[9] = $seq; # rewrite the sequence

				# get quality value if wished
				if($quality) {
					$qual = <$fh_fastq>; # description line
					$qual = <$fh_fastq>; # quality line
					chomp $qual;

					if($reverse == 1) {
						$qual = reverse($qual);
						$reverse = 0;
					}
					if($cstart != 0) {
						$qual = substr($qual, $cstart);
						$cstart = 0;
					}
					if($cend != 0) {
						$qual = substr($qual, 0, $cend);
						$cend = 0;
					}
					$tmp[10] = $qual; # rewrite the quality
				}

				# update use tag of that fastq sequence
				$usedID = $data[1];
				$mappingStatus[$usedID][1] = 1;
			}
			else {
				print "not found ID: $id\n";
				$notFound++;
			}
			# print the new line
			print $fh_out join("\t", @tmp)."\n";
		}
		else {
			print "invalid split count at line \n$row\n";
			close($filenameOUT);
			unlink $filenameOUT;
			exitCode("EXIT_MISFORMATED_INPUT");
		}
		$i++;

		# status update
		if($i % 1000000 == 0) {
			print "processed $i sam entries ($notFound sequenes are missing)...\n";
		}	
	}	
}
print "finished parsing of SAM file ($i entries, $notFound sequenes are missing)!\n";

if($notFound > 0) {
	close($filenameOUT);
	unlink $filenameOUT;
	print "$notFound sequenes are missing!\n";
	exitCode("EXIT_MISFORMATED_INPUT");
}

# write unmapped reads to disk
if($filenameUnmapped ne "-1") {
	print "started writing unmapped reads to FASTQ file...\n";
	open($fh_unmapped, '>:encoding(UTF-8)', $filenameUnmapped) or do { print "Could not open output file '$filenameUnmapped'!\n"; exitCode("EXIT_WRITING_FAILED"); };
	my $ii = 0;
	my $iw = 0;
	my @data;

	# check all entries
	for($i = 0; $i < scalar(@mappingStatus); $i++) {
		@data = @{$mappingStatus[$i]};

		# check if sequence was not aligned
		if($data[1] == 0) {
			# select correct file handle
			$fh_fastq = $fileHandles[$data[2]];
			# jump to correct position
			seek($fh_fastq, $data[0], 0);

			# read four lines
			for($ii = 0; $ii < 4; $ii++) {
				$row = <$fh_fastq>;  
				chomp $row;  
				print $fh_unmapped $row."\n"; 
			}
			$iw++; 
		}
		# status update
		if($i % 1000000 == 0) {
			print "processed $i fastq entries ($iw are unmapped)...\n";
		} 
	}
	print "finished writing of unmapped FASTQ file ($iw entries)!\n";
} 

# close file handles
close($fh_sam);
foreach $fh_fastq (@fileHandles) {
	close($fh_fastq);
}
print "Finished!\n";

if($filenameUnmapped ne "-1") {
	close($fh_unmapped);
}

# test if output file is really there
if (-e $filenameOUT) {
	if($filenameUnmapped ne "-1") {
		if(-e $filenameUnmapped) {
			writeReturnValues($filenameOUT, $filenameUnmapped);
			exitCode("EXIT_OK");
		}
		else {
			exitCode("EXIT_FAILED");
		}
	}
	writeReturnValues($filenameOUT, "");
	exitCode("EXIT_OK");
} 
else {
	exitCode("EXIT_FAILED");
}

sub revCom {
	my $dna = $_[0];
	$dna = reverse($dna);
	$dna =~ tr/ACGTacgt/TGCAtgca/;
	return $dna;
}

sub writeReturnValues {
	if($returnFilePath ne "") {
		my $samFile = $_[0];
		my $unmappedFile = $_[1];
		open(my $fh_outVars, '>:encoding(UTF-8)', $returnFilePath) or do { print "Could not open return param file '$returnFilePath'!\n"; exitCode("EXIT_RETURN_PARAMETER_MISSING"); };
		print $fh_outVars "SAMFileWithSequences\t$samFile\n";
		print $fh_outVars "UnmappedReadFile\t$unmappedFile\n";
		print $fh_outVars "?EOF!\n";
		$fh_outVars->autoflush;
		close($fh_outVars); 
	}
}

sub getReadIDWithSuffix {
	$id = $_[0];
	$flag = $_[1];

	# strip ID at blank
	my @tmp = split(" ", $id);
	$id = $tmp[0];
		
	# add /1 or /2 if paired read
	if($flag & 64) {
		$id="$id/1"
	}
	elsif($flag & 128) {
		$id="$id/2"
	}
	else {
		$id="$id/1"
	}
	return $id;
}
