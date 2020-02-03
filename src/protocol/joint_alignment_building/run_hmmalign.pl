#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin;
use lib $FindBin::Bin . '/../../lib';
use Localpaths qw(load_localpaths);

my %files;


my $local_paths 	= load_localpaths();
my $project_root 	= $local_paths->{root_dir};
my $src_dir 		= $project_root . $local_paths->{src_dir};
my $cdhit_dir 		= $project_root . $local_paths->{cdhit};
my $hmmalign_path	= $project_root . $local_paths->{hmmer}  . "/binaries/hmmalign";

my $input_file1				= "";
my $input_file2				= "";
my $stockholm_output_file1	= "";
my $stockholm_output_file2	= "";
my $fasta_output_file1		= "";
my $fasta_output_file2		= "";
my $num_cpus 				= 1;
my $hmm1 					= "";
my $hmm2 					= "";

GetOptions (
	"input_file1=s" 	=> \$input_file1,
	"input_file2=s" 	=> \$input_file2,
	"sto1=s"		 	=> \$stockholm_output_file1,
	"sto2=s"		 	=> \$stockholm_output_file2,
	"fasta1=s" 			=> \$fasta_output_file1,
	"fasta2=s" 			=> \$fasta_output_file2,
	"hmm1=s" 			=> \$hmm1,
	"hmm2=s" 			=> \$hmm2,
);


# TODO: Check parameters one by one and define usage subroutine
if(!$input_file1 || !$input_file2 || !$stockholm_output_file1 || !$stockholm_output_file2 || !$fasta_output_file1 
	|| !$fasta_output_file2 || !$hmm1 || !$hmm2){

	die "ERROR while running run_hmmalign.pl. A mandatory parameter have not been defined\n";
}

# Run hmmalign command for the two input files, one for each pfam domain
run_hmmalign($input_file1, $input_file2, $stockholm_output_file1, $stockholm_output_file2, $hmmalign_path, $hmm1, $hmm2);

# Reformat files form Pfam stockholm format to FASTA
reformat($stockholm_output_file1, $stockholm_output_file2, $fasta_output_file1, $fasta_output_file2);


exit;



# This subroutine call hmmalign on alingment files (they start with PF) in input_dir and run hmmalign for the first and second domain
# num_cpus indicates the number of hmmalign processes that are run on parallel
sub run_hmmalign{

	my $input_file1				= shift;
	my $input_file2 			= shift;
	my $stockholm_output_file1 	= shift;
	my $stockholm_output_file2 	= shift;
	my $hmmalign_path 			= shift;
	my $hmm1 					= shift;
	my $hmm2 					= shift;
	
#	print("Aligning file $input_file1 againts $hmm1 Pfam HMM profile\n");
	print("Aligning sequences againts $hmm1 HMM profile\n");
	system("$hmmalign_path --allcol --trim --outformat Pfam -o $stockholm_output_file1 $hmm1 $input_file1");
#	print("Aligning file $input_file2 againts $hmm2 Pfam HMM profile\n");
	print("Aligning sequences againts $hmm2 HMM profile\n");
	system("$hmmalign_path --allcol --trim --outformat Pfam -o $stockholm_output_file2 $hmm2 $input_file2");
}


# This subroutine change the format of any alignment in input dir from Pfam stockholm to FASTA (insertions are removed, only HMM matches states are saved)
sub reformat{
	
	my $stockholm_output_file1 	= shift;
	my $stockholm_output_file2 	= shift;
	my $fasta_output_file1 		= shift; 
	my $fasta_output_file2 		= shift;

	# TODO: Change this piece of code for something much cleaner. Process by perl internally
	print "Reformating alignments files from Stockholm to FASTA\n";
	open(FH, $stockholm_output_file1) || die "Error while opening file $stockholm_output_file1";
	open(FH_out, ">$fasta_output_file1") || die "Error while opening file $fasta_output_file1";
#	my $str = 'grep -v -P "^(#|\n|\s|\/)" ' . "$stockholm_output_file1 " . ' | awk \'/^.+/ {print ">"$1"\n"$2}\' ' . " > $fasta_output_file1";
#	`$str`;
	while(my $line = <FH>){
		if($line =~ /^#/){
			next;
		}
		elsif($line =~ /^\//){
			next;
		}
		elsif($line =~ /^$/){
			next;
		}
		else{
			my @fields = split(/\s+/, $line);
			print FH_out ">$fields[0]\n$fields[1]\n";
		}
	}
	close FH;
	close FH_out;

#	print "Reformating file $stockholm_output_file2\n";
	open(FH, $stockholm_output_file2) || die "Error while opening file $stockholm_output_file2";
	open(FH_out, ">$fasta_output_file2") || die "Error while opening file $fasta_output_file2";
	while(my $line = <FH>){
		if($line =~ /^#/){
			next;
		}
		elsif($line =~ /^\//){
			next;
		}
		elsif($line =~ /^$/){
			next;
		}
		else{
			my @fields = split(/\s+/, $line);
			print FH_out ">$fields[0]\n$fields[1]\n";
		}
	}
	close FH;
	close FH_out;
}


