#!/bin/bash

# Evaluate the kmer the distribution of unique single copy kmers in minimap alignments compared to the reads.
# We hope to show that:
	# 1) showing that reads with exactly have only one alignment have unique kmers (this approach is not losing anything)
	# 2) for reads that have multiple good alignments, minimap’s primary choice is not always the same as the one that our 
	#    approach would pick (this approach is better for picking correct alignments)
	# 3) Reads that minimap drops contain unique kmers that can be used to map (this approach allows us to align more 
	#    reads accurately)


# minimap SAM output, including the read itself as the last column in the file for each row/alignment
mult_align_file=$1
ref_unique_db=$2

# Get all the first best alignments and make a new meryl database
# sort by name? Primary is denoted with the flag 0x900, supplementary with 0x800? chimeric?


samtools view -f0x0 $mult_align_file | samtools bam2fq
samtools view -f0x100 $mult_align_file > supplementary.bam


# Retrieve the reads from the bam files to create the databases
samtools view primary.bam | awk $1;$3 (mapping score) \n $7 (the read sequence) > primary_reads.fasta
samtools view supplementary.bam | awk $1;$3 (mapping score) \n $7 (the read sequence) > supplementary_reads.fasta

# Make the databases

meryl count k=15 primary_reads.fasta output primary.meryl
meryl count k=15 supplementary_reads.fasta output supplementary.meryl

# intersect with the reference sck database

meryl intersect primary.meryl $ref_unique_db output primary_intersect.meryl
meryl histogram primary_intersect.meryl > primary_intersect.histo

meryl intersect supplementary.meryl $ref_unique_db supplementary_intersect.meryl
meryl histogram supplementary_intersect.meryl > supplementary_intersect.histo