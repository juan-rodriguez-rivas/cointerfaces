#!/usr/bin/env perl
use strict;
use Getopt::Long;
use FindBin;
use lib $FindBin::Bin . '/../../../lib';
use Localpaths qw(load_localpaths);

my $local_paths 			= load_localpaths();
my $tmp_dir 				= $local_paths->{temporary_directory};

my $input_file;
my $output_file;

GetOptions (
	"input_file=s" 	=> \$input_file,
	"output_file=s" => \$output_file,
);

if(!$input_file){
	print "ERROR: Input file has to be defined";
	usage();
}
if(!$output_file){
	print "ERROR: Output file has to be defined";
	usage();
}

# Check input file exists and has non zero size
if(! -s $input_file){
	die "ERROR: Input file non exists or is empty";
}

if($input_file eq $output_file){
	die "ERROR: Input and output files can not be the same file";	
}

# Remove those pairs that contains at least one sequence that appear more than once in the alignment
remove_duplications($input_file, $output_file);

exit;

# Remove those pairs that contains at least one sequence that appear more than once in the alignment
sub remove_duplications{

	my $input_file 	= shift;
	my $output_file = shift;
	
	my $count = 0;
	my %remove_gene;
	my %boundaries;

	print "Removing ambiguous sequences\n";

	open(FH, $input_file) || die "Error while opening file $input_file";
	
	while(my $line = <FH>){
		
		if($line =~ /^>/){
			
			my @genes = split(/\t/, $line);
			
			my @gene1 = split(/_/, $genes[0]);
			my @gene2 = split(/_/, $genes[1]);
			
			my $gene1 = substr($gene1[0], 1);
			my $gene2 = substr($gene2[0], 1);
			
			$remove_gene{dom1}{$gene1}{count}++;
			$remove_gene{dom2}{$gene2}{count}++;
			
			if(!defined $remove_gene{dom1}{$gene1}{line}){
				$remove_gene{dom1}{$gene1}{line} = [$count];
			}
			else{
				push(@{$remove_gene{dom1}{$gene1}{line}}, $count);
			}
			
			if(!defined $remove_gene{dom2}{$gene2}{line}){
				$remove_gene{dom2}{$gene2}{line} = [$count];
			}
			else{
				push(@{$remove_gene{dom2}{$gene2}{line}}, $count);
			}				
		}
		
		$count++;
	}
	
	my %remove_lines;
	foreach my $dom (keys %remove_gene){
		foreach my $gen (keys %{$remove_gene{$dom}}){
			# Remove from the alignment those sequences that appear more than once on the alignment of its domain
			if($remove_gene{$dom}{$gen}{count} > 1){
				foreach my $count (@{$remove_gene{$dom}{$gen}{line}}){
					$remove_lines{$count} 	= 1;
					$remove_lines{$count+1} = 1;						
				}
			}
			
			# TODO: Remove from the alignment those pairs that appear in both orders
			# Remove from the alignment those pairs that appear in both orders (due to paralogous that have a significant 
			# score in both pfam HMM profiles)
			# We have to check whether one of them is a substring of the another
			# EDIT: Not necessary now, paraloguos will overlap so overlapping filters will deal with them
		}
	}
	
	
	$count = 0;
	seek(FH, 0, 0);
	
	open(FH_OUT, ">$output_file") || die "Error while opening for writing output file $output_file";
	while(my $line = <FH>){
		
		if($remove_lines{$count}){
			$count++;
			next;
		}
		
		print FH_OUT "$line";
		
		$count++;
	}
	close FH_OUT;
}



