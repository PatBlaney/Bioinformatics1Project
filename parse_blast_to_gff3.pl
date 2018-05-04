#!/usr/bin/perl
use warnings;
use strict;
use feature 'say';

## parse_blast_to_gff3.pl
## Written by C Roesel 2013, modified by Verena 2016.

## Modified by Patrick Blaney 2018 as per Module 12 assignment

# This program runs a BLAST then parse the output into a GFF3 formatted file.

# Run BLAST of the CRISPR sequences stored in fasta against the publically accessible database
# DrosophilaAllChroms
my($perfect, $off_target) = blast_kmers( "unique_crisprs.fasta", "DrosophilaAllChroms" );

# Print perfect BLAST hits to GFF3 file in proper format
print_gff3_file("crisprs_with_ot.gff3", $perfect, $off_target);

close GFF3FILE;

sub blast_kmers {
	my ( $input_fasta, $db ) = @_;

	# Put BLAST command and all params into an array. This could be created as single
	# string, but an array makes it easier to see individual parameters.
	my @blast_command = (
		'blastn', # run nuc/nuc BLAST program
		'-task blastn', # change task from 'megablast' to 'blastn'
		'-db ' . $db, # database (subject sequence)
		'-query ' . $input_fasta, # query sequences
		'-outfmt 6' # output in tabular format
	);

	# Print BLAST command for debugging purposes.
	say "@blast_command";

	# Execute BLAST command and retrieve output in BLAST filehandle.
	open( BLAST, "@blast_command |" );
	
	# Process BLAST output line-by-line. Ignore lines that start with '#' (comments).
	while ( my $blast_output = <BLAST> ) {
		chomp $blast_output;
		
		# Initialize variables to hold the hash references that are returned at end of
		# subroutine
		my ($perfect_hits, $off_target_hits) =
			process_blast_output( $blast_output ) if $blast_output !~ /^#/;

		# Test if hashes carried over
		#print keys(%$off_target_hits);
		#print keys(%$perfect_hits);
		#print values(%$off_target_hits);
		#print values(%$perfect_hits);
	}
}

sub process_blast_output {
	my ( $blast_output_line ) = @_;

	# Because output is in tabular format, line can be split by tabs. Each column value
	# is then assigned to meaningful variables.
	my (
		$query_id, $chrom, $pident, $length, $mismatches,
		$gaps, $q_start, $q_end, $start, $end
	) = split( "\t", $blast_output_line );
	
	# Initialize hash to store perfect BLAST hits
	my %perfect_BLAST_hits = ();
	
	# Initialize a count that will help keep track of occurances of the perfect BLAST hits
	my $hit_count = 1;
	
	# Check if the identity is 100% and the sequence has a length of 21, if so 
	if ($pident == 100 && $length == 21) {
		
		# Check if perfect BLAST has already been stored, if not then store it as
		# a key with the current index as the value. 
		if (not defined $perfect_BLAST_hits{$blast_output_line}) {

			$perfect_BLAST_hits{$blast_output_line} = $hit_count;
		}
		else {
			$perfect_BLAST_hits{$blast_output_line} += 1;
		}
	}
	# Initialize hash to store count of each off target
	my %off_targets = ();
	
	# Initialize a the start count of each off target
	my $off_target_count = 1;

	# Keep count of all off target BLAST hits, they have an identify less than 100% and
	# a mismatch count of 3 or fewer
	if ($pident < 100 && ($mismatches + (21 - $length)) <= 3) {
		if (not defined $off_targets{$blast_output_line}) {
			
			# Store entry into hash with initial count
			$off_targets{$blast_output_line} = $off_target_count;
		}
		else {
			# Since entry already exists in hash, increment the count by 1
			$off_targets{$blast_output_line} += 1;
		}
	}

	# Return the perfect hits and off target hits as hash references
	return \%perfect_BLAST_hits, \%off_targets;
}

sub print_gff3_file {
	my ($filename_to_write, $BLAST_hash_ref, $off_target_hash_ref) = @_;
	
	# Open .gff3 file to write perfect BLAST hits in GFF3 format to
	open(GFF3FILE, ">", $filename_to_write);
	
	# Use foreach loop to move through all BLAST entries which were stored as keys
	foreach my $BLAST_key (keys %$BLAST_hash_ref) {
		 
		# Test
		#print $BLAST_key, "\n";

	        # Because BLAST input is in tabular format, line can be split by tabs.
        	# Each column value is then assigned to meaningful variables.
        	my (	$query_id, $chrom, $pident,
			$length, $mismatches,
			$gaps, $q_start, $q_end,
			$start, $end )
					 = split( "\t", $BLAST_key );

        	# Change local list separator to tabs.
        	local $, = "\t";
  
        	# Set default strand as forward strand.
        	my $strand = '+';
        	my $gff_start = $start;
        	my $gff_end = $end;
 
	       	# Overwrite to reverse strand if start position is greater than end position.
      		if ( $start > $end ) {
              		$strand = '-';
              		$gff_start = $end;
              		$gff_end = $start;
      		}
 	
		# ..... 
		#my $note_update = $$off_target_hash_ref{$query_id};
	
      		my @new_line = (
              		$chrom, ".", 'OLIGO', $gff_start, $gff_end, ".", $strand, ".",
              		"Name=$query_id;Note= note_update off-target matches"
     		);
		
		print GFF3FILE  @new_line, "\n";
	}	
}
