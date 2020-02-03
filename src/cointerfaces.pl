#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd;
use FindBin;
use lib $FindBin::Bin . '/lib';
use Localpaths qw(load_localpaths);
use File::Basename;
use List::Util qw(shuffle);


my $local_paths 			= load_localpaths();
my $project_root 			= $local_paths->{root_dir};
my $tmp_dir 				= $project_root . $local_paths->{temp_dir};
my $src_dir 				= $project_root . $local_paths->{src_dir};

my $base_pairs_operon_distance	= 300;
my $strategy_uniques 			= 1;

my $query_pfams_file;
my $hmm1;
my $hmm2;

my $output_root;
my $output_ali1;
my $output_ali2;
my $output_scores;
my $working_dir = cwd();


GetOptions (
	"hmm1=s" 				=> \$hmm1,
	"hmm2=s" 				=> \$hmm2,
	"genomic_distance=f" 	=> \$base_pairs_operon_distance,
	"output_root=s"			=> \$output_root,
	"tmp_dir=s"          	=> \$tmp_dir,
);


## Test
#$hmm1 = "PF00005.22.hmm";
#$hmm2 = "PF00528.17.hmm";

check_input($hmm1, $hmm2);

# If not provided, the output_root is derived from the working dir and the filenames of the input HMM profiles
if(!defined $output_root){
	$output_root = $working_dir . "/$hmm1-$hmm2";  
}
$output_ali1 = $output_root . ".fasta1";
$output_ali2 = $output_root . ".fasta2";
$output_scores = $output_root . ".zscores.inter";
`mkdir -p $tmp_dir`;

# Using the HMM profiles, search homologus sequences in the genomes and pair those in genomic proximity (or unique in the genome if this strategy is enable)
# As output generates a join alignment (split in two fasta files)
build_joint_alignments($hmm1, $hmm2, $output_ali1, $output_ali2, $strategy_uniques, $base_pairs_operon_distance);

# From the joint alignment, computes the DCA coevolutionary model
compute_dca_model($hmm1, $hmm2, $output_ali1, $output_ali2, $output_scores);


exit;

# This function generates a joint alingments compose by two alingments in FASTA whose sequences have been paired by genomic distance (and uniqness in genomes if this option is enabled)
sub build_joint_alignments{
	my $hmm1 = shift;
	my $hmm2 = shift;
	my $output_ali1 = shift;
	my $output_ali2 = shift;
	my $strategy_uniques = shift;
	my $base_pairs_operon_distance = shift;
	
	# Search the genomes for homologous sequences using the HMM profiles
	my $hmmsearches_dir = $tmp_dir. "/hmmsearches";
	system("$src_dir/protocol/joint_alignment_building/search_pfams_in_genomes.pl -hmm1 $hmm1 -hmm2 $hmm2 -output_dir $hmmsearches_dir");
	
	# Find those pairs of domain sequences that are likely interacting
	my $pairing_dir = $tmp_dir. "/pairing";
	system("$src_dir/protocol/joint_alignment_building/make_pairs.gene_neighbouring.uniques.pl -hmm1 $hmm1 -hmm2 $hmm2 -hmmsearches $hmmsearches_dir -output_dir $pairing_dir");
	
	# Check if it has been found something. Otherwise, it finishes
	my $root_pairing_file = "$pairing_dir/$hmm1-$hmm2";
	my $pairing_file1 = "$root_pairing_file.fasta1";
	my $pairing_file2 = "$root_pairing_file.fasta2";
	if(-z $pairing_file1 && -z $pairing_file2){
		print "There is not any pair of interacting pair of sequence domains using the HMM profiles: $hmm1 and $hmm2\n";
		exit;
	}
	
	# Align the sequences against their Pfam profiles and apply the quality filters
	my $joint_alignment_dir = $tmp_dir. "/joint_alignment";
	system("$src_dir/protocol/joint_alignment_building/align_and_filter.pl -input_file1 $pairing_file1 -input_file2 $pairing_file2 -output_dir $joint_alignment_dir -output_file1 $output_ali1 -output_file2 $output_ali2 -hmm1 $hmm1 -hmm2 $hmm2");
	
	# Check if the output_file exist and are non-empty (for some cases, no joint sequences will be found or will lose due to the filters. No error output, it is normal)
	if(! -s $output_ali1 || ! -s $output_ali2){
		print "There are not any sequences after applying quality filters in case $hmm1\t$hmm2\n"; 
		exit;
	}
	
	print "Joint alingment built\n";
}


sub check_input{
	my $hmm1 = shift;
	my $hmm2 = shift;
	
	if(!defined $hmm1){
		print "\nERROR: A pair of first HMM profile (-hmm1 hmm_file1 option) has to be provide\n";
		usage();
	}
	
	if(!defined $hmm2){
		print "\nERROR: A pair of first HMM profile (-hmm2 hmm_file2 option) has to be provide\n";
		usage();
	}

	if($hmm1 eq $hmm2){
		print "\nERROR: hmm1 and hmm2 should be different\n";
		usage();
	}	
}


# It receives a joint alignment (two fasta files), prepare and compute mpl program (fortran implementation of plmDCA), i.e. the DCA coevolutionary model
sub compute_dca_model{
	
	my $hmm1 			= shift;
	my $hmm2 			= shift;
	my $input_file1 	= shift;
	my $input_file2 	= shift;
	my $output_file 	= shift;

	my $temp_dca_dir = $tmp_dir. "/dca_model";
	`mkdir -p $temp_dca_dir`;
	
	my $hmm1_basename = basename($hmm1);
	my $hmm2_basename = basename($hmm1);
	my $temp_file = "$temp_dca_dir/$hmm1-$hmm2.mpl";

	# Preprocess the joint alignment
	system("$src_dir/protocol/dca_model/tools/make_mpl_input.py -1 $input_file1 -2 $input_file2 -o $temp_file");
	#if(! -s $out_file){print STDERR "No output_file found while preprocessing joint alignment to compute the coevolutionary model\n"; return}
	
	# Compute the coevolutionary model for joint alignment
	system("$src_dir/protocol/dca_model/mpl/mpl -i $temp_file -l 0.01 -g");
	#if(! -s $out_file){print STDERR "No output_file found while computing the coevolutionary model\n"; return}
	
	# Renumber the output files to meet the HMMs numbering and apply APC correction
	system("$src_dir/protocol/dca_model/tools/renumber_scores.py -s $temp_file.scores -d $temp_file.dim -o $temp_file");
	
	# Apply MAD standarization
	system("$src_dir/protocol/dca_model/tools/z-scores.py -i $temp_file.scores.inter -o $temp_file.zscores.inter -f 3 > $temp_file.zscores.inter");

	# Copy to final file
	`cp $temp_file.zscores.inter $output_file`;
	
	
	print "DCA model computed\n";
	print "Cointerfaces protocol ended\n";
}





sub usage{
	
	print "\nUsage:\n";
	print "./cointerfaces.pl -hmm1 HMM_file1 -hmm2 HMM_file2\n";
	print "Example:\n";
	print "./run_protocol.pl -hmm1 PF00177.hmm -hmm2 PF00380.hmm\n";
	print "\nDescripion:\n";
	print "Build a joint multiple sequence of likely interacting domain sequences and computes a coevolutionary model from the alignment. Highly coeolving pairs of residues are usually contacts between the domains in the 3D structure. If the output root is not specified, they are generated from the filenames of the HMM profiles. Three files are generate, two fasta files codifying the joint alignment a file with output APC zscores for each pair of position of the two HMM profiles\n";
	print "\nParameters:\n";
	print "-hmm1 pfam_id -> pfam identifier of the first putatively interacting.\n";
	print "-hmm2 pfam_id -> pfam identifier of the second putatively interacting.\n";
	print "[-output_root]       -> Optional. Root of the output files\n";
	print "[-genomic_distance]  -> Optional. Maximum genomic distance (in base pairs) to consider two domain in two different genes as interacting. By default, 300 base pairs is used. Overrides the value specified in the localpaths file\n";
	print "[-temp_dir]          -> Optional. Directory where the intermediate files will be stored. Overrides the value specified in the localpaths file\n";
	

	exit; 
}
