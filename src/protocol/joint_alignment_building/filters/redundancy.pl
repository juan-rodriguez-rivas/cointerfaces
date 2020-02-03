#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
use FindBin;
use lib $FindBin::Bin . '/../../../lib';
use Localpaths qw(load_localpaths);


my $local_paths 	= load_localpaths();
my $project_root 	= $local_paths->{root_dir};
my $cdhit_dir 		= $project_root . $local_paths->{cdhit};

# Mandatory parameters
my $input_file;
my $output_file;
my $output_dir;

# Optional parameters
my $redundancy = 0.8;


GetOptions (
	"input_file=s" 		=> \$input_file,
	"output_file=s" 	=> \$output_file,
	"output_dir=s" 		=> \$output_dir,
	"redundancy=f" 		=> \$redundancy,
);

# Check input file exists and has non zero size
if(! -s $input_file){
	die "ERROR: Input file non exists or is empty";
}

if($input_file eq $output_file){
	die "ERROR: Input and output files can not be the same file";	
}


find_no_redundant($input_file, $output_file, $redundancy);


exit;



sub find_no_redundant{
	
	my $input_file	= shift;
	my $output_file	= shift;
	my $redundancy	= shift;
	
	print "Removing sequences with high redundancy\n";

	`$cdhit_dir/cd-hit -c $redundancy -i $input_file -o $output_file`;
}

