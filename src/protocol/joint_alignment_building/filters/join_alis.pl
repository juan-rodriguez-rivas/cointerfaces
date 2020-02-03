#!/usr/bin/perl -w
use strict;

# Simple script that given two fasta files, concatenate both (with a tab in the middle)

my $input_file1 = $ARGV[0];
my $input_file2 = $ARGV[1];
my $output_file = $ARGV[2];

open(FH1, $input_file1) || die "Error while opening file $input_file1";
open(FH2, $input_file2) || die "Error while opening file $input_file2";
open(FH_OUT, ">", $output_file) || die "Error while opening file $output_file";
my @lines1 = <FH1>;
my @lines2 = <FH2>;
chomp @lines1;
chomp @lines2;
my $num_lines = scalar @lines1;
for (my $i=0; $i < $num_lines; $i++){
	print FH_OUT "$lines1[$i]\t$lines2[$i]\n";
}
close FH1;
close FH2;
close FH_OUT;


exit;