package Localpaths;
use strict;
use warnings;
use base 'Exporter';
use File::Basename;

our @EXPORT = qw/ load_localpaths /;

sub load_localpaths{
	
	my $local_paths_file = dirname(__FILE__) . "/../localpaths";
	
	my %paths;
	open(FH, $local_paths_file) || die "Error while opening file localpaths";
	while(my $line = <FH>){
		
		if($line !~ /^\w/){
			next;
		}
		
		chomp $line;
		
		my @fields = split(/\s*=\s*/, $line);
		
		my $data = $fields[0];
		my $path = $fields[1]; 
		$paths{$data} = $path;
	}
	close FH;
	
	return \%paths;
}


1;