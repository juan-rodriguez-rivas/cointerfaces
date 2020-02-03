#!/usr/bin/perl -w
#use List::MoreUtils qw/ uniq /;
use strict;
use Getopt::Long;
use FindBin;
use lib $FindBin::Bin . '/../../lib';
use Localpaths qw(load_localpaths);
use Cwd;
use File::Basename;

# Load config
my $local_paths 	= load_localpaths();
my $project_root 	= $local_paths->{root_dir};
my $genomes_list 	= $project_root . $local_paths->{list_of_genomes};
my $genomes_fasta 	= $project_root . $local_paths->{ensembl_bacteria_fasta};
my $genomes_gtf 	= $project_root . $local_paths->{ensembl_bacteria_gtf};
my $genomes_info	= $project_root . $local_paths->{genomes_info};
my $hmmalign_bin 	= $project_root . $local_paths->{hmmer} . "/binaries/hmmalign";
my $base_pairs_operon_distance	= $local_paths->{operon_distance};
my $strategy_uniques 			= $local_paths->{strategy_uniques};
my $output_dir 					= cwd();

# Parameters for internal functionalities
# If $select_strains=1, only one strain per species if used. Otherwise, all available genomes are used
my $select_strains				= 0;
my $count_failed_sequences		= 0;
my $num_pairs 					= 0;

# Get input parameters
my $hmm1;
my $hmm2;
my $hmmsearches_dir;
GetOptions (
	"hmmsearches_dir=s"		=> \$hmmsearches_dir,
	"genomic_distance=f" 	=> \$base_pairs_operon_distance,
	"output_dir=s" 			=> \$output_dir,
	"hmm1=s" 				=> \$hmm1,
	"hmm2=s" 				=> \$hmm2,
	"strategy_uniques=f"	=> \$strategy_uniques,
);

## Test
#$hmm1 = "/home/jrodrig5/projects/cointerfaces/test/PF00005.22.hmm";
#$hmm2 = "/home/jrodrig5/projects/cointerfaces/test/PF00528.17.hmm";
#$hmmsearches_dir = "/home/jrodrig5/projects/cointerfaces/tmp/hmmsearches";
#$output_dir = "/home/jrodrig5/projects/cointerfaces/tmp/joint_alignment";


if(!defined $hmm1 ||  !defined $hmm2){
	print STDERR "\nERROR: A pair of HMM profiles (-hmm1 hmm_id1 and -hmm2 hmm2_id options) have to be provide\n";
	usage();
}

# Variables for mapping genomes and taxids and their hierarchy
my %get_taxid;
my %get_specie_name;
my %get_parent;
my %get_sons;
my @species;
my %species;

load_species($genomes_list);
load_taxid($genomes_list, $genomes_info);


my $count_domains_instances = 0;
my %interactions;
my %search_loci;
my %genes_loci_and_sequence;


# Recover HMM files basename, hmmsearches are stored according to them
my $hmm1_basename = basename($hmm1);
my $hmm2_basename = basename($hmm2);

# Load hmmsearches results. Every gene found is marked to be included in pairing process. 
# Starting and ending position of the hits are stored
print "Loading hmmsearches\n"; 
load_hmmsearch_results($hmm1_basename, $hmm2_basename, $genomes_list, $hmmsearches_dir);

# Load the genomic coordinates of the genes that were found using hmmsearch
print "Loading genomic coordinates of genes\n";
load_loci_and_sequence($genomes_gtf, $genomes_fasta);

my %select_strain;
my %pairs;
print "Pairing sequences\n";
make_pairing($hmm1_basename, $hmm2_basename);
print_alis($hmm1_basename, $hmm2_basename, $output_dir);
print "Number of paired sequences = " . (scalar @{$pairs{$hmm1_basename}{$hmm2_basename}}) . "\n";



exit;



sub load_loci_and_sequence{
	
	my $gtf_files = shift;
	my $pep_files = shift;
	
	my %gene_id_genome;
	my %sequences_genome;
	
	my $num_genomes = scalar keys %get_taxid;
	
	# Store temporary the location of those genes that we have found through the hmmsearch saved in %genes_loci_and_sequence
	foreach my $genome (sort keys %get_taxid){
		
		undef %gene_id_genome;
		
		# Find the genomic region annotation for the genes
		my $file_in = `ls $gtf_files/$genome*`;
		chomp $file_in;
		if(!$file_in || ! -e $file_in){print "File $file_in not found. Genome $genome\n"; next}
		open(FH, $file_in) || die "Error while opening file $file_in";

		my $start;
		my $end;
		my $pot_gene_id;
		my $count 	= 1;
		my $contig 	= '';
		while(my $line = <FH>){
			
			my @fields = split(/\s+/, $line);
			
			$pot_gene_id = $fields[11];
			
			$pot_gene_id =~ s/[";]//g;
			
			if(!defined $genes_loci_and_sequence{$genome}{$pot_gene_id}){
				next;
			}
			
			# Check if the gene appear more than once in the genome. In this case, it is discarded (untrusted annotation)
			if(defined $gene_id_genome{$pot_gene_id}){
				print STDERR "Warning: Gene $pot_gene_id found more than once in genome $genome, line $count\n";
				next;
			}
			
			$contig			= $fields[0];
			$start 			= $fields[3];
			$end 			= $fields[4];
			
			$gene_id_genome{$pot_gene_id} = 1;
			$genes_loci_and_sequence{$genome}{$pot_gene_id}{contig}	= $contig;
			$genes_loci_and_sequence{$genome}{$pot_gene_id}{start} 	= $start;
			$genes_loci_and_sequence{$genome}{$pot_gene_id}{end} 	= $end;
			$genes_loci_and_sequence{$genome}{$pot_gene_id}{count}++;
			
			$count++;
		}
		close FH;
		
	
		# Find the FASTA sequence for each gene	
		undef %sequences_genome;
		
		$file_in = `ls $pep_files/$genome*`;
		chomp $file_in;
		if(!$file_in || ! -e $file_in){print "Fasta genomic file $file_in not found. Genome $genome\n"; next}
		open(FH, $file_in) || die "Error while opening file $file_in";

		my $gene_id = '';
		my $new_gene_id = '';
		my $sequence;
		my $jump = 0;
		while(my $line = <FH>){
			
			if($line =~ /^>(\w+)\s/){
				
				$new_gene_id = $1;
				if($gene_id && $sequence){
					$genes_loci_and_sequence{$genome}{$gene_id}{sequence} = $sequence;
				}
				$sequence = '';
				
				if(!defined $genes_loci_and_sequence{$genome}{$new_gene_id}){
					$gene_id = '';
					$jump = 1;
				}
				else{
					$gene_id = $new_gene_id;
					$jump = 0;
				}
				
			}
			else{
				if($jump){
					next;
				}
				
				chomp $line;
				$sequence .= $line;
			}
		}
		
		# For the last gene in the genome
		if($gene_id && $sequence){
			$genes_loci_and_sequence{$genome}{$gene_id}{sequence} = $sequence;
		}		
		

		close FH;
		
	}
	
	# Check things are ok or delete the gene otherwise
	foreach my $genome (sort keys %get_taxid){
		
		foreach my $gene_id (keys %{$genes_loci_and_sequence{$genome}}){
			
			my $start 		= $genes_loci_and_sequence{$genome}{$gene_id}{start};
			my $sequence 	= $genes_loci_and_sequence{$genome}{$gene_id}{sequence};
			
			if( !$start || !$sequence){
					
				delete $genes_loci_and_sequence{$genome}{$gene_id};
				
				print STDERR "Warning: No chromosome loci found for gene $gene_id in genome $genome\n";
			}
		}
	}
}




sub print_alis{
	
	my $hmm1 = shift;
	my $hmm2 = shift;
	my $output_dir = shift;
	`mkdir -p $output_dir`;
	
			
	my %selected_strains;
	
	my $hmm_file1 = "$output_dir/$hmm1-$hmm2.fasta1";
	my $hmm_file2 = "$output_dir/$hmm1-$hmm2.fasta2";
	
	
	open(FH_hmm1, ">$hmm_file1");
	open(FH_hmm2, ">$hmm_file2");
	
	my $count = 1;
	foreach my $par (@{$pairs{$hmm1}{$hmm2}}){
				
		my $taxid 	= $par->[0];
		
		if($select_strains){
			
			my $parent = $get_parent{$taxid};
			
			if(!defined $parent){
				$selected_strains{$taxid} = $taxid;
			}
			
			if($selected_strains{$parent} != $taxid){
				next;
			}
		}
		
		my $genome 		= $get_specie_name{$taxid};
		my $gene_id1 	= $par->[1];
		my $gene_id2 	= $par->[2];
		my $overlap 	= $par->[3];
		my $sequence1 	= $genes_loci_and_sequence{$genome}{$gene_id1}{sequence};
		my $sequence2 	= $genes_loci_and_sequence{$genome}{$gene_id2}{sequence};
		
		if(!defined $sequence1 || !defined $sequence2){
		    next;
		}
		
		# Sequences with '*' are problematic, they are discarded
		if($sequence1 =~ /\*/){
			next;
		}
		
		if($sequence2 =~ /\*/){
			next;
		}
		
		my $start1			= $interactions{$hmm1}{$genome}{$gene_id1}{domain_start};
        my $end1			= $interactions{$hmm1}{$genome}{$gene_id1}{domain_end};
		my $difference1 	= $end1 - $start1 + 1;
				
		my $start2 			= $interactions{$hmm2}{$genome}{$gene_id2}{domain_start};
		my $end2			= $interactions{$hmm2}{$genome}{$gene_id2}{domain_end};
		my $difference2 	= $end2 - $start2 + 1;


		if($overlap){
			# If there some overlap between two domains of the same protein, we cut one of the domain sequence to avoid repeat the same residues
			# Otherwise, correlations will be perfect between those positions
			# If it positive we cut the first residues of the second domain. Otherwise, at the beginning of the first domain
			if($overlap > 0){
				$start2 = $start2 + $overlap;
			}
			elsif($overlap < 0){
				$start1 = $start1 - $overlap;
			}
		}
		
		my $sequence1_domain = substr($sequence1, $start1-1, $difference1);
		my $sequence2_domain = substr($sequence2, $start2-1, $difference2);
		
		my $length1 = length($sequence1_domain);
		my $length2 = length($sequence2_domain);
		
		if($length1 == 0 || $length2 == 0){
			$count_failed_sequences++;
			next;
		}
		
		my $gene_start1 = $genes_loci_and_sequence{$genome}{$gene_id1}{start};
		my $gene_start2 = $genes_loci_and_sequence{$genome}{$gene_id2}{start};
		my $gene_end1 	= $genes_loci_and_sequence{$genome}{$gene_id1}{end};
		my $gene_end2 	= $genes_loci_and_sequence{$genome}{$gene_id2}{end};
		
		my $genomic_distance = 0;
		if($gene_id1 ne $gene_id2){
			
			my $max_start = $gene_start1;
			if($gene_start2 > $gene_start1){
				$max_start = $gene_start2;
			}
			
			my $min_end = $gene_end2;
			if($gene_end1 < $gene_end2){
				$min_end = $gene_end1;
			}
			
			$genomic_distance = $max_start - $min_end;
		}						
		
		if($length1 > 0 && $length2 > 0){
			print FH_hmm1 ">$gene_id1\_$taxid\/$start1-$end1\n$sequence1_domain\n";
			print FH_hmm2 ">$gene_id2\_$taxid\/$start2-$end2\n$sequence2_domain\n";
#			print FH_hmm2 ">$gene_id2\_$taxid\/$start2-$end2\t$genomic_distance\n$sequence2_domain\n";
		}
	}
	close FH_hmm1;
	close FH_hmm2;
	
	# If in some case, there is no pairs, we want to create an empty output files to show that there no pairs in this case
	
	if(! -e $hmm_file1 && ! -e $hmm_file2){
		open(FH_hmm1, ">$hmm_file1");
		open(FH_hmm2, ">$hmm_file2");
		close FH_hmm1;
		close FH_hmm2;		
	}
}



# This subroutine finds the pairs of sequence from both domains that meet the conditions of likely interactions:
# - They are in the same protein but are not overlapping or overlap only a small part
# - If they are in two different proteins:
#   -- They are uniques in the whole genome (both domains appear just once in the whole genome). Only if option uniques is activated.
#   -- They are in the same operon, i.e. their genomic distance is below a given threshold.
sub make_pairing{

	my $hmm1 = shift;
	my $hmm2 = shift;

#	my $num_cases 		= $num_pairs * (scalar keys %get_taxid);
			
	foreach my $genome (keys %{$interactions{$hmm1}}){
		
		foreach my $gene_id1 (keys %{$interactions{$hmm1}{$genome}}){
			
			if(	!defined $genes_loci_and_sequence{$genome}{$gene_id1}){
				next;
			}
			
			my $locus_start1	= $genes_loci_and_sequence{$genome}{$gene_id1}{start};
			my $locus_end1		= $genes_loci_and_sequence{$genome}{$gene_id1}{end};
			
			foreach my $gene_id2 (keys %{$interactions{$hmm2}{$genome}}){
				
				# If strategy uniquess is included:
				# We check if both domains appear just once in the genome, then they are paired
				if($strategy_uniques){
					
					if($interactions{$hmm1}{$genome}{count} == 1 && $interactions{$hmm2}{$genome}{count} == 1 && $gene_id1 ne $gene_id2){
						my @array = ($get_taxid{$genome}, $gene_id1, $gene_id2, 0);
						if(defined $pairs{$hmm1}{$hmm2}){
							push(@{$pairs{$hmm1}{$hmm2}}, @array);
						}
						else{
							$pairs{$hmm1}{$hmm2} = [\@array];
						}
						
						next;
					}
				}
				
				if(!defined $genes_loci_and_sequence{$genome}{$gene_id2}){
					next;
				}
						
				my $contig1			= $genes_loci_and_sequence{$genome}{$gene_id1}{contig};
				my $contig2			= $genes_loci_and_sequence{$genome}{$gene_id2}{contig};
				
				if($contig1 ne $contig2){
					next;
				}
						
				my $locus_start2	= $genes_loci_and_sequence{$genome}{$gene_id2}{start};
				my $locus_end2		= $genes_loci_and_sequence{$genome}{$gene_id2}{end};
				
				my $domain_start1   = $interactions{$hmm1}{$genome}{$gene_id1}{domain_start};
				my $domain_end1     = $interactions{$hmm1}{$genome}{$gene_id1}{domain_end};

				my $domain_start2   = $interactions{$hmm2}{$genome}{$gene_id2}{domain_start};
				my $domain_end2     = $interactions{$hmm2}{$genome}{$gene_id2}{domain_end};


				# If both domains belong to the same gene, it checks that there not too much overlap between them 
				# We allow a maximum overlap of the 10% of the smaller domain (HMM profile length) to a maximum of 20 residues 
				# Si es el mismo gen, comprobar que no hay overlap entre los dominios
				if($gene_id1 eq $gene_id2){
							
					# First it checks when there is no overlap
                    if( $domain_end1 < $domain_start2 ||
                    	$domain_end2 < $domain_start1){
										
						my @array = ($get_taxid{$genome}, $gene_id1, $gene_id2, 0);
						if(defined $pairs{$hmm1}{$hmm2}){
							push(@{$pairs{$hmm1}{$hmm2}}, \@array);
						}
						else{
							$pairs{$hmm1}{$hmm2} = [\@array];
						}
						next;									
					}
					
							
					# If they do not overlap, it checks that there not too much overlap between them
					my $tam_dom1		= $interactions{$hmm1}{$genome}{$gene_id1}{length_hmm}*0.1;
					my $tam_dom2		= $interactions{$hmm2}{$genome}{$gene_id2}{length_hmm}*0.1;
					
					my $overlap = 0;
					if($tam_dom1 > $tam_dom2){
						$overlap = $tam_dom2;
					}
					else{
						$overlap = $tam_dom1;
					}
					
					if($overlap > 20){$overlap = 20}
							
                    if( $domain_end1 < $domain_start2 + $overlap ||
                    	$domain_end2 < $domain_start1 + $overlap){
									
						# Compute the exact overlap
						my $exact_overlap = 0;
						if($domain_start2 > $domain_start1){
							# First domain appears first in the genome
							if($domain_end1 > $domain_end2){
								# It should not occur, the second domain is embebbed in the first domain
								next;										
							}
							$exact_overlap = $domain_end1 - $domain_start2 + 1; 
						}
						else{
							if($domain_end2 > $domain_end1){
								# It should not occur, the second domain is embebbed in the first domain
								next;										
							}
							# El domain 2 es el primero, negativo indica que el segundo dominio es el primero en el genoma
							$exact_overlap = -($domain_end2 - $domain_start1) - 1;
						}
						
						my @array = ($get_taxid{$genome}, $gene_id1, $gene_id2, $exact_overlap);
						if(defined $pairs{$hmm1}{$hmm2}){
							push(@{$pairs{$hmm1}{$hmm2}}, \@array);
						}
						else{
							$pairs{$hmm1}{$hmm2} = [\@array];
						}
						next;									
					}
					else{
						# Too much overlap, it is discarded
						next;
					}
				}
				else{
					
					# It the two domains not belong to the same gene, it checks that the genomic distance between the two genes
					# is smaller than the given "operon distance" threshold
					if(abs($locus_start1 - $locus_end2) < $base_pairs_operon_distance ||
					   abs($locus_start2 - $locus_end1) < $base_pairs_operon_distance){
						
						
						my @array = ($get_taxid{$genome}, $gene_id1, $gene_id2, 0);
						if(defined $pairs{$hmm1}{$hmm2}){
							push(@{$pairs{$hmm1}{$hmm2}}, \@array);
						}
						else{
							$pairs{$hmm1}{$hmm2} = [\@array];
						}
					}
				}
			}
		}
	}
}


# Load hmmsearch results, all the sequence hits will be annotate in $genes_loci_and_sequence including where the hit starts and end 
sub load_hmmsearch_results{
	
	my $hmm1 			= shift;
	my $hmm2 			= shift;
	my $genomes_file 	= shift;
	my $dir 			= shift;
	
	
	my @hmms = ($hmm1, $hmm2);
	
	
	# Load hmmsearches in each genome for each HMM. Looking first at genomes allow us to take advantage of cached hard disk lectures 
	foreach my $genome (keys %get_taxid){	
		foreach my $hmm (@hmms){
			
			my $file_in = "$dir/$genome/$hmm";
			
			if(!(-e $file_in)){
				next;
			}
			
			my @in;
			open(FH, $file_in) || die "Error while opening file $file_in";
			while(my $line = <FH>){
				if($line =~ /^#/){
					next;
				}
				chomp $line;
				push(@in, $line);
			}
			close FH;
			
			foreach my $domain_instance (@in){
				
				my @fields = split(/\s+/, $domain_instance);
				
				# If the domain appears more than once on the gene, they are removed as we can not be sure which pair 
				# will interact or not
				my $count_domain = $fields[10];
				if($count_domain > 1){
					$interactions{$hmm}{$genome}{count} += $count_domain;
					next;
				}
				
				$count_domains_instances++;
				$interactions{$hmm}{$genome}{count}++;
				my $new_domain = \%{$interactions{$hmm}{$genome}};
				
				# Mark this gene to be annotated afterwards in %genes_loci_and_sequence. Store starting and ending position of the match in %interactions
				my $gene_id = $fields[0];
				$new_domain->{$gene_id}{length_hmm} 	= $fields[5];
				$new_domain->{$gene_id}{domain_start} 	= $fields[17];
				$new_domain->{$gene_id}{domain_end} 	= $fields[18];
				if(!defined $genes_loci_and_sequence{$genome}{$gene_id}){
					$genes_loci_and_sequence{$genome}{$gene_id} = {};
				}
				else{
					
				}
				
				# Include this gene for the search of possible pairs
				my @array = [$gene_id, $hmm];
				if(defined $search_loci{$genome}){
					push(@{$search_loci{$genome}}, @array);
				}
				else{
					$search_loci{$genome} = [\@array];
				}
			}
		}
	}
}



sub load_sons{
	
	my @appear;
	my @not_appear;
	
	foreach my $taxid (keys %get_specie_name){
		
		if(defined $get_parent{$taxid}){
			
			push(@appear, $taxid);
			
			if(defined $get_sons{$get_parent{$taxid}}){
				my @array = @{$get_sons{$get_parent{$taxid}}};
				push(@array, $taxid);
				$get_sons{$get_parent{$taxid}} = @array;
			}
			else{
				$get_sons{$get_parent{$taxid}} = [$taxid];
			}
		}
		else{
			push(@not_appear, $taxid);
		}
	}
}



sub load_parents{
	
	my %dev;
	my $file_in = shift;
	my $taxonomy_file;
	open($taxonomy_file, $file_in) || die "Error while opening file taxonomy";
	
	my $id 		= '';
	my $parent	= '';
	my $count 	= 0;
	while(my $line = <$taxonomy_file>){
		
		if($line =~ /^ID\s+:\s+(\d+)/){
			$id = $1;
		}
	
		if($line =~ /^PARENT ID\s+:\s+(\d+)/){
		
			if(!defined $id){
				print STDERR "Error in taxonomy file in line $count: Undefined ID\nline: $line\n";
			}
			else{
				$parent = $1;
	
				if(defined $get_parent{$id}){
					print STDERR "Previously defined parent in taxonomy: line $count\n";
				}
				else{
					if(defined $get_specie_name{$id}){
						$get_parent{$id} = $parent;
					}
				}
			}
		}
		
		if($line =~ /^\/\//){
			$id 	= '';
			$parent = '';
		}
		
		$count++;
	}
	
	close $taxonomy_file;
}


sub load_species{
	
	my $file_in = shift;
	
	open(FH, $file_in) || die "Error while opening file $file_in";
	
	while(my $line = <FH>){
		
		if($line =~ /^\n/){
			next;
		}
		
		chomp $line;
		$species{$line} = 1;
	}
	
	close FH;
}


sub load_taxid{
	
	my $genomes_list = shift;
	my $genomes_info = shift;

	
	my %all_taxids;
	open(FH, $genomes_info) || die "Error while opening file $genomes_info";
	while(my $line = <FH>){
		my @fields = split("\t", $line);
		my $assembly = $fields[5];
		my $taxid = $fields[3];
		
		$all_taxids{$assembly} = $taxid;
	}
	close FH;
	
	open(FH, $genomes_list) || die "Error while opening file $genomes_list";
	
	while(my $genome_name = <FH>){
		chomp $genome_name;
		
		my $assembly;
		my $taxid;
		if($genome_name =~ /.*(GCA_\d+\.\d*)/){
			$assembly = $1;
			
			if($all_taxids{$assembly}){
				$get_taxid{$genome_name} = $all_taxids{$assembly};
				$get_specie_name{$all_taxids{$assembly}} = $genome_name;
			}
			else{
				print "Unrecognized assembly in genome info\n";
			}		
		}
		else{
			print "Unrecognized assembly in genome list\n";
		}
		
	}
	close FH;
}



sub usage{
	
	print "\nUsage:\n";
	print "./make_pairs.gene_neighbouring.pl -hmm1 hmm_id -hmm2 hmm_id\n";
	print "Example:\n";
	print "./make_pairs.gene_neighbouring.pl -hmm1 PF10417.4 -hmm2 PF00578.16\n";
	print "\nDescripion:\n";
	print "From the results of a profiles finds all the likely interacting pairs of domains sequences across prokaryote genomes. Results are stored in two FASTA files, one for each hmm domain. The two files have the same number of lines and the pairs of sequences with the same number of line in both files are very likely interacting protein domains.\n";
	print "\nParameters:\n";
	print "-hmm1 hmm_id -> hmm identifier of the first putatively interacting domain whose pairs have to be calculated.\n";
	print "-hmm2 hmm_id -> hmm identifier of the second putatively interacting domain whose pairs have to be calculated.\n";
	print "\tPF10417.4\tPF00578.16\n\tPF10417.4\tPF02195.13\n\tPF00244.15\tPF00583.19\n";
	print "[-output_dir] -> Output directory where the output files will be stored. By default, the working directory (the one where the program is called)\n\n";

	exit; 
}





