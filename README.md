# Bioinformatics1Project
Final homework/project for Bioinformatics Methods 1

The goal of the program is to take a FASTA formatted file, perform a BLAST over the sequences located in the file, parse the BLAST output lines for perfect matches, write those matches to a txt file in gff3 format

Feedback from submission: 
"Rather than declare two new hashes each time process_blast_output() is called (which is every time a BLAST output line is read in) and adding key-value pairs into them, consider declaring the hashes in blast_kmers() instead, passing them into the subroutine call, and adding key-value pairs into those variables"

"Keys used for hashes are too unique -- there are no duplicate BLAST hits in the BLAST output. here are, however, duplicate crispr IDs, which might be better used as the key. That way, you will be able to ask, "What is the perfect BLAST hit (hash 1's value) for this crispr ID (key) and how many off-targets (hash 2's value) were found for it (key)?"

"Rather than keep count of the number of occurrences of perfect BLAST hits, consider saving the converted BLAST hit as the value for %perfect_BLAST_hits. That way, you don't have to re-process the line again in print_gff3_file(), a.k.a. re-split the line to convert it to gff3 format."
