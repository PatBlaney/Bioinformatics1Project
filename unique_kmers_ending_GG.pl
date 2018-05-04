#!/usr/bin/perl
use strict;
use warnings;

# Store all sequences found in the file passed to the subroutine within an array
my $sequence = load_sequence("/scratch/Drosophila/dmel-all-chromosome-r6.02.fasta");

# Sanity check to see that sequences were properly stored and are retrievible
#print $$sequence[1869], "\n";


# Using the sequences from the file that were stored in the array, populate one reference hash with
# the full kmer sequence as the keys and the last 12 positions of that respective kmer sequence as the
# value. Populate another reference hash with the last 12 positions of each kmer sequence as the keys
# and the count associated with the occurance of that sequence 
my ($kmers_count, $full_kmer_sequences) = count_kmers($sequence, 21);

# Sanity check to see if hashes were populated properly and are accessible 
#print keys(%$kmers_count), "\n";
#print $$kmers_count{"TATTCGACCTGG"}, "\n";

#print keys(%$full_kmer_sequences), "\n";
#print $$full_kmer_sequences{"TATTCGACCTGG"}, "\n";

my $output_filename = 'unique_crisprs.fasta';
print_kmers_to_fasta($kmers_count, $full_kmer_sequences, $output_filename);

sub load_sequence {
	my ($seq_file) = @_;

	# Create an empty array to populate with the list of sequences from the file
	my @seq_set = ();

	open (SEQ, "<", $seq_file) or die "Cannot open file: $!\n";


	# Initialize a count to index the sequences into the array of all sequences in the file
	# Start is -1 to populate the 0 index of the array
	my $count = -1;
	
	while (<SEQ>) {
		chomp;
		
		# If line is header, increment count to move the index since this means a new
		# sequence is going to be read
		if ($_ =~ /^>/) { 
			$count += 1;
		}
		
		# Add the read in lines to the array of sequences
		else {
			$seq_set[$count] .= $_;
		}
			
	}
	
	# For testing
	#print $count, "\n";
	
	close SEQ;

	# Return a reference array of all sequences within the file
	return \@seq_set;
}

sub count_kmers {
	my ($seq_array, $kmer_length) = @_;
	
	# Create an empty hash to populate with the last 12 positions of each unique k-mer sequence
	# and count associated with each sequence
	my %kmer_count = ();

	# Create a second empty hash to populate with the the 12-end sequence and its associated full
	# length sequence
	my %full_kmer_seqs = (); 
	
	# Utilize a foreach loop to iterate through each sequence in the array of all seqences
	# from the file, nest for loop to extract k-mers
	foreach my $seq (@$seq_array) {
		for (my $seq_index = 0; $seq_index <= length($seq) - $kmer_length; $seq_index++) {
		
			# Extract each k-mer of length 21 from the sequence provided as the total
			# sequence
			my $kmer_sequence = substr($seq, $seq_index, $kmer_length);
			
			# For testing
			#print $kmer_sequence, "\n";

			# Sanity check to see if sequence is only 21 nucleotides long and the
			# sequence ends with two guanines. Also used to separately captrue the
			# only the last 12 positions of the sequence 
			if ($kmer_sequence =~ /([ATGC]{9}([ATGC]{10}GG))/) {
				
				# Capture full kmer sequence that matches 5'-NGG-3' pattern
				my $full_kmer_sequence = $1;

				# Capture last 12 positions of kmer sequence
				my $last_12_nucleotides = $2;
				
				# Given each new kmer sequence an initial count if not already found in
				# the k-mer count hash
				my $count = 1;
			
				# First check to see if the last 12 sequence has already been stored
				# into the count hash
				if (not defined $kmer_count{$last_12_nucleotides}) {
					
					# Store the separately captrued last 12 positions of the kmer
					# sequence as a key in the kmer count hash with associated
					# count as the value
					$kmer_count{$last_12_nucleotides} = $count;

					# Store the full kmer sequence into the alternate hash as a
					# value with the associated sequence of the last 12 positions
					# as the key
					$full_kmer_seqs{$last_12_nucleotides} = $full_kmer_sequence;
				}
				else {	
					# If the full kmer sequence has already been found and stored
					# as a key in the alternate hash, then simply increment the
					# count associated with the last 12 position sequence
					$kmer_count{$last_12_nucleotides} += 1;
				}
			}
		}
	}
	# Return both hash references 
	return \%kmer_count, \%full_kmer_seqs;
}

sub print_kmers_to_fasta {
	my ($crispr_count_ref, $full_kmer_seq_ref, $file_to_write) = @_;
	open (KMERS, ">", $file_to_write);
	
	# Initialize a count to give each crispr entry and index within the header
	my $count = 1;

	# Cycle through keys of hash containing crispr sequences to find all unique sequences then
	# write the full kmer sequence of these unique crisprs to a .txt file
	foreach my $crispr_key (keys %$crispr_count_ref) {
		
		# Using the count associated with each crispr sequence, identify if the sequence is
		# unique. If it is, print the full kmer sequence to the .txt file
		if ($$crispr_count_ref{$crispr_key} == 1) {
			
			# Creat the entry within the .txt file with the header that includes the index
			# of each entry. If the crispr key is unique, use it as the key for the 
			# reference hash to output the full kmer sequence to the .txt file
			print KMERS ">crispr_$count\n$$full_kmer_seq_ref{$crispr_key}\n";
			
			# Increment the count after each entry
			$count += 1;
		}
	}
	close KMERS;
}		
