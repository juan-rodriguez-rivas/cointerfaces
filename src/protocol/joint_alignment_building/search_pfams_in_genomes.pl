#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;
use lib $FindBin::Bin . '/../../lib';
use Localpaths qw(load_localpaths);
use Getopt::Long;
use File::Basename;

my $hmm1;
my $hmm2;
my $output_dir;

GetOptions (
	"hmm1=s"		=> \$hmm1,
	"hmm2=s"		=> \$hmm2,
	"output_dir=s"	=> \$output_dir,
);

## Test
#$hmm1 = "/home/jrodrig5/projects/cointerfaces/test/PF00005.22.hmm";
#$hmm2 = "/home/jrodrig5/projects/cointerfaces/test/PF00528.17.hmm";

my $local_paths 	= load_localpaths();
my $project_root 	= $local_paths->{root_dir};
my $genomes_list 	= $project_root . $local_paths->{list_of_genomes};
my $genomes_fasta 	= $project_root . $local_paths->{ensembl_bacteria_fasta};
my $hmmsearch_bin 	= $project_root . $local_paths->{hmmer} . "/binaries/hmmsearch";

`mkdir -p $output_dir`;

#TODO: Check files and dirs exists
run_hmmsearch($hmm1, $hmm2, $genomes_list, $genomes_fasta, $output_dir);

exit;


sub run_hmmsearch{
	
	my $hmm1 = shift;
	my $hmm2 = shift;
	my $genomes_list = shift;
	my $genomes_fasta = shift;
	my $output_dir = shift;
	
	my @hmms = ($hmm1,$hmm2);
	my @genomes = `cat $genomes_list`;
	chomp @genomes;
	
	for my $hmm (@hmms){ 
		print "Searching for homologous sequences using the HMM profile: $hmm\n";
	
	foreach my $genome (@genomes){

		my $genome_path = `ls $genomes_fasta/$genome*pep*`;
		my $output_genome_path = "$output_dir/$genome";
		
		my $basename_hmm = basename($hmm);
		`mkdir -p $output_genome_path`;
    	my $output_file_temp 	= "$output_dir/$genome/$basename_hmm\_temp";
    	my $output_file_final 	= "$output_dir/$genome/$basename_hmm";
    		if(! -e $output_file_final){
    			my $cmd = "$hmmsearch_bin --noali --domtblout $output_file_temp --cut_ga $hmm $genome_path";
	    		`$hmmsearch_bin --noali --domtblout $output_file_temp --cut_ga $hmm $genome_path`;
	    		`grep -v "^#" $output_file_temp > $output_file_final`;
   				`rm $output_file_temp`;
    		}
		}
	}
}
