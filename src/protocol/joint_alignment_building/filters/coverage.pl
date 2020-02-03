#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
use FindBin;
use lib $FindBin::Bin . '/../../../lib';
use Localpaths qw(load_localpaths);

# Mandatory parameters
my $input_file;
my $output_file;
my $output_dir;

# Optional parameters
my $coverage = 0.8;


GetOptions (
	"input_file=s" 		=> \$input_file,
	"output_dir=s" 		=> \$output_dir,
	"output_file=s" 	=> \$output_file,
	"coverage=f" 		=> \$coverage,
);

# Check input file exists and has non zero size
if(! -s $input_file){
	die "ERROR: Input file non exists or is empty";
}

if($input_file eq $output_file){
	die "ERROR: Input and output files can not be the same file";	
}


check_coverage($input_file, $output_file, $output_dir, $coverage);


exit;



sub check_coverage{
	
	my $input_file	= shift;
	my $output_file	= shift;
	my $tmp_dir 	= shift; 
	my $coverage	= shift;
	
	print "Removing sequences with low coverage\n";

	open(FH, ">$output_file") || die "Error while opening file $output_file $!";
	
	my %coverage;
	my %remove_lines;

	my $input_filename = basename($input_file);

	my $input_file1 = $tmp_dir . "/" . $input_filename . ".cov1";
	my $input_file2 = $tmp_dir . "/" . $input_filename . ".cov2";
	
#	my $output_filename = basename($output_file);
#
#	my $output_file1 = $tmp_dir . "/" . $output_filename . ".cov1";
#	my $output_file2 = $tmp_dir . "/" . $output_filename . ".cov2";
	
	open(FH, $input_file) 			|| die "Error while opening file $input_file: $!";
	open(FH1, ">$input_file1") 		|| die "Error while opening file $input_file1: $!";
	open(FH2, ">$input_file2") 		|| die "Error while opening file $input_file2: $!";
	my @lines = <FH>;
	foreach my $line (@lines){
		my @fields = split(/\t/, $line);
		print FH1 $fields[0] . "\n";
		print FH2 $fields[1];
	}
	close FH;
	close FH1;
	close FH2;
		
	# For HMM 1		
	check_domain_coverage($input_file1, \%remove_lines);
	
	# For HMM 2
	check_domain_coverage($input_file2, \%remove_lines);
	
	open(FH, $input_file) 			|| die "Error while opening file $input_file: $!";
	open(FH1, "$input_file1") 		|| die "Error while opening file $input_file1: $!";
	open(FH2, "$input_file2") 		|| die "Error while opening file $input_file2: $!";
	open(FH_OUT, ">$output_file") 	|| die "Error while opening file $output_file: $!";
		
	my $count = 0;
	# All files have the same number of lines
	while(my $line = <FH>){
		my $line1 = <FH1>;
		my $line2 = <FH2>;
		
		if(!$remove_lines{$count}){
			
			if($line1 =~ /^>/){
				print FH_OUT "$line";
			}
			else{
				chomp $line1;
				print FH_OUT "$line1\t$line2";
			}
		}
		
		$count++;
	}
	
	close FH1;
	close FH2;
	close FH_OUT;
	
	`rm $input_file1`;
	`rm $input_file2`;
	
	close FH;
}




sub check_domain_coverage{
	
	my $file			= shift;
	my $remove_lines 	= shift;
	
	my %coverage;
	
	
	my $max = 0;
	
	open(FH, $file) || die "Error while opening file $file: $!";
	
	my $header;
	
	# The length of the profile is the number of mayuscules plus the number of hyphens in any sequence (we just take the first)
	my $line = <FH>; # Fasta header, discard
	
	$line =~ /^>(.+)/;	
	$header = $1;
	$line = <FH>; # This is the sequence

	# Find the profile
	my $result = 0;
	$result++ while($line =~ m/\p{Uppercase}/g);
	$result++ while($line =~ m/-/g);
	
	$max = $result;
	
	# We save the number of hits for the first domain
	$result = 0;
	$result++ while($line =~ m/\p{Uppercase}/g);
	$coverage{$header} = $result;

	# Count the number os mayuscules of each sequence to discard them or not later on	
	while(my $line = <FH>){
		
		if($line =~ /^>(.+)/){
			$header = $1;
		}
		else{
			
			chomp $line;
			my $result = 0;
			# Only count mayuscules, no gaps, this what we want to remove
			$result++ while($line =~ m/\p{Uppercase}/g);
			$coverage{$header} = $result;
		}
	}
	
	my $threshold = $max * $coverage;
	
	# Restart from the beggining and discard sequences with low coverage (many gaps)
	seek(FH, 0, 0);
	
	my $file_name = substr($file, 0, -8);
	
	my $count = 0;
	$header = '';
	while(my $line = <FH>){
		
		if($line =~ /^>(.+)/){
			
			$header = $1;
			
			if(!defined $coverage{$header}){
				print "Warning: No coverage found for sequence $coverage{$header})\n";
				$count++;
				next;
			}
			
			if($coverage{$header} < $threshold){
				$remove_lines->{$count} 	= 1;
				$remove_lines->{$count+1} 	= 1;
			}
		}

		$count++;
	}		
	
	close FH;
}

