#!/usr/bin/perl -w
use strict;
use File::Basename;
use Cwd;
use Getopt::Long;
use FindBin;
use lib $FindBin::Bin . '/../../lib';
use Localpaths qw(load_localpaths);

my $local_paths 	= load_localpaths();
my $project_root 	= $local_paths->{root_dir};
my $src_dir 		= $project_root . $local_paths->{src_dir};
my $hmmer_dir 		= $project_root . $local_paths->{hmmer};
my $hhsuite_dir		= $project_root . $local_paths->{hhsuite};


# Mandatory parameters
my $input_file1;
my $input_file2;
my $output_file1;
my $output_file2;
my $output_dir;
my $hmm1;
my $hmm2;

# Optional parameters
my $redundancy 				= 0.8;
my $coverage 				= 0.8;


GetOptions (
	"input_file1=s" 		=> \$input_file1,
	"input_file2=s" 		=> \$input_file2,
	"output_dir=s" 			=> \$output_dir,
	"output_file1=s" 		=> \$output_file1,
	"output_file2=s" 		=> \$output_file2,
	"hmm1=s" 				=> \$hmm1,
	"hmm2=s" 				=> \$hmm2,
	"coverage=f" 			=> \$coverage,
	"redundancy=f" 			=> \$redundancy,
);


## Test
#$hmm1 = "/home/jrodrig5/projects/cointerfaces/test/PF00005.22.hmm";
#$hmm2 = "/home/jrodrig5/projects/cointerfaces/test/PF00528.17.hmm";
#$input_file1 = "/home/jrodrig5/projects/cointerfaces/tmp/joint_alignment/PF00005.22.hmm-PF00528.17.hmm.fasta1";
#$input_file2 = "/home/jrodrig5/projects/cointerfaces/tmp/joint_alignment/PF00005.22.hmm-PF00528.17.hmm.fasta2";
#$output_dir = "/home/jrodrig5/projects/cointerfaces/tmp/joint_alignment";

if(!$input_file1 || !$input_file2){
	print "\nERROR: Two input files and two output files have to be specified\n";
	usage();
}


# All four input files have to be different. Check if that is the case
my @files;
push(@files, $input_file1, $input_file2);
check_files(\@files);


my $hmm1_basename = basename($hmm1);
my $hmm2_basename = basename($hmm2);
my $filename1 = "$output_dir/$hmm1_basename";
my $filename2 = "$output_dir/$hmm2_basename";
my $filename_both = "$output_dir/$hmm1_basename-$hmm2_basename";

#my $output_file1 	=  "$output_dir/$hmm1_basename.ali1";
#my $output_file2 	=  "$output_dir/$hmm1_basename.ali2";

# Needed intermediate files
my $file_aligned_sto1 		= "$filename1.sto";
my $file_aligned_sto2 		= "$filename2.sto";
my $file_aligned_fasta1 	= "$filename1.fasta";
my $file_aligned_fasta2 	= "$filename2.fasta";
my $file_just_profile1 		= "$filename1.profile";
my $file_just_profile2 		= "$filename2.profile";
my $file_joint		 		= "$filename_both.joint";
my $file_no_amb 			= "$filename_both.no_amb";
my $file_cov 				= "$filename_both.cov";
my $file_red 				= "$filename_both.red";


if(! -s $input_file1){
	print "Error in align_and_filter.pl: Input file $input_file1 does not exists or is empty\n";
}
if(! -s $input_file2){
	print "Error in align_and_filter.pl: Input file $input_file2 does not exists or is empty\n";
}

# Align the sequence against their Pfam profiles and apply all the filter step by step checking that output files exists and are not empty
system("$src_dir/protocol/joint_alignment_building/run_hmmalign.pl -input_file1 $input_file1 -input_file2 $input_file2 -sto1 $file_aligned_sto1 -sto2 $file_aligned_sto2 -fasta1 $file_aligned_fasta1 -fasta2 $file_aligned_fasta2 -hmm1 $hmm1 -hmm2 $hmm2");
if(! -s $file_aligned_fasta1 || ! -s $file_aligned_fasta1){exit;}
system("$src_dir/protocol/joint_alignment_building/filters/just_profile.bash $file_aligned_fasta1 $file_aligned_fasta2 $file_just_profile1 $file_just_profile2");
if(! -s $file_just_profile1  || ! -s $file_just_profile2){exit;}
system("$src_dir/protocol/joint_alignment_building/filters/join_alis.pl $file_just_profile1 $file_just_profile2 $file_joint");
if(! -s $file_joint){exit;}
system("$src_dir/protocol/joint_alignment_building/filters/remove_ambiguities.pl -input_file $file_joint -output_file $file_no_amb");
if(! -s $file_no_amb){exit;}
system("$src_dir/protocol/joint_alignment_building/filters/coverage.pl -input_file $file_no_amb -output_dir $output_dir -output_file $file_cov -coverage $coverage");
if(! -s $file_cov){exit;}
system("$src_dir/protocol/joint_alignment_building/filters/redundancy.pl -input_file $file_cov -output_file $file_red -redundancy $redundancy");
if(! -s $file_red){exit;}
system("$src_dir/protocol/joint_alignment_building/filters/split_alis.bash $file_red $output_file1 $output_file2");
if(! -s $output_file1 || ! -s $output_file2){exit;}


exit;

sub usage{
	
	print "\nUsage:\n";
	print "./align_and_filter -hmm1 hmm1_file -hmm1 hmm2_file -input_file1 file1 -input_file2 file2 -output_dir dir -output_file1 output_file1 -output_file2 output_file1\n";
	print "\nDescription:\n";
	print "From two files with pairs of interacting sequences of two HMMs, align both files to their corresponding pfam profiles and apply filters to remove ambiguity, coverage and redundancy. The four files porvided have to be different\n";
	print "\nParameters:\n";
	print "-hmm1 -> File of the first HMM profile.\n";
	print "-hmm1 -> File of the second HMM profile.\n";
	print "-input_file1 -> File with the sequences of the first HMM domain.\n";
	print "-input_file1 -> File with the sequences of the second HMM domain.\n";
	print "-output_dir  -> Directory where intermediate results will be stored";
	print "-output_file1 -> File where the resulting sequences of the first domain will be stored.\n";
	print "-output_file2 -> File where the resulting sequences of the second domain will be stored.\n";
	print "-redundancy -> Value (from 0 up to 1) of the redundancy filter corresponding to the proportion of identical amino acids (sequence identity). Overrides the valued set in localpaths\n";
	print "-coverage   -> Value (from 0 up to 1) of the coverage filter corresponding to the proportion of amino acids (no gaps) covering the profile (e.g. 0.8 mean 80% of the profile covered, a maximum of 20% of gaps per domain sequence). Overrides the valued set in localpaths\n";
	
	exit;
}



sub check_files{
	my $files = shift;
	
	for(my $i = 0; $i < scalar @$files; $i++){
		for(my $j = $i+1; $j < scalar @$files; $j++){
			if($files->[$i] eq $files->[$j]){
				print STDERR "\nERROR: All four input and output files have to be different. Repeated file: $files[$i]\n";
				usage();
				exit;
			}
		}
	}
}
