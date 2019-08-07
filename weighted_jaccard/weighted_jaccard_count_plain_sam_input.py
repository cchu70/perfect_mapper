# Script to count the shared and unique kmers between two sequences
# copy of weighted_jaccard_count.py, but I use a normal sam file

from  weighted_jaccard_func import getKmers, parseKmerFile, parseFasta, counts, weightJaccard, parseCigar
import sys
import os
import subprocess
import time
import threading
import numpy

# python ../../scripts/mashmap_postfilter/weighted_jaccard/weighted_jaccard_count_plain_sam_input.py GAGE_A.sim_reads.fasta error_0.0/AAF.err_0.0_A.v_1.fasta error_0.0/AAF_minimap2.N50_r3k.split.err_0.0_A.v_1.aligned_A.sam k_file 21 GAGE_A GAGE_A 0.0


def parseSam(sam_str):
	data = sam_str.split()

	read_name = data[0]

	ref_name = data[2]

	ref_start = int(data[3])
	cigar = data[5]

	length, read_start, read_end = parseCigar(cigar)

	if length < 0:
		sys.stderr.write("Unable to parse cigar string %s. Original string was: \n%s\n" % (cigar, sam_str))	
		

	ref_end = ref_start + length

	return read_name, length, ref_name, ref_start, ref_end, read_start, read_end


k_size = 21


def main():

	read_fasta = sys.argv[1]				# Reads used to align to the ref fasta file, producing the sam file
	ref_fasta = sys.argv[2]					# Ref fasta file that the reads above mapped to, producing the sam file. This particular ref fasta introduced some err_rate to one of the partitions (which_error)
	sam_file = sys.argv[3]					# The aformentioned sam files
	k_file = sys.argv[4]					# the kmer list of exisiting kmers, formatted kmer_seq<tab>frequency, and all frequencies must be > 0 
	k_size = int(sys.argv[5])				# The size of the kmers used

	which_reads_aligned = sys.argv[6]		# Name of the reads used to align (ex. "GAGE_A")
	which_error = sys.argv[7]				# Name of the partition in the reference that we simualted error (ex. "GAGE_A")
	err_rate = float(sys.argv[8])			# The error rate introduced for book keeping purposes
	
	# only care if these are the same
	# if which_reads_aligned != which_error:
	# 	# Only consider the events where the reads we are mapping originates from the region we introduced error
	# 	return
	# #####


	unique_kmer_weight = float(sys.argv[9])
	non_unique_kmer_weight = float(sys.argv[10])


	# Get read sequences
	sys.stderr.write("Parsing Read fasta: %s\n" % read_fasta)
	read_records = parseFasta(read_fasta) # Dictionary of read names and it's corresponding sequence

	sys.stderr.write("Parsing Ref fasta: %s\n" % ref_fasta)
	ref_record = parseFasta(ref_fasta)

	sys.stderr.write("Parsing kmer file: %s\n" % k_file)
	kmer_table = parseKmerFile(k_file)
	sys.stderr.write("Finished counts kmers\n")


	# Initialize Loop
	curr_read_str = None
	curr_read_name = None


	max_score = (0, "")

	correct_count = 0
	incorrect_count = 0
	unmapped_count = 0

	for line in open(sam_file, "r"):

		if "@" not in line:		# skipping header

			read_name, length, ref_name, ref_start, ref_end, read_start, read_end =  parseSam(line.strip())

			if length < 0:
				sys.stderr.write("Original sam file was: \n%s\n" % (sam_file))
				unmapped_count += 1
			else:

				# Check which read (current or next) this alignment corresponds to 
				if (curr_read_str):
					if (read_name != curr_read_name):
						# evaluate the curr read performance
						# print("%s" % max_score[1])

						if max_score[1]:
							# the current read actually had an alignment
							#sys.stderr.write("Max score on %s, aligned %s\n" % ( max_score[1], which_reads_aligned))
							if max_score[1] == which_reads_aligned: # GAGE_B
								#sys.stderr.write("Correct! %s = %s" % (max_score[1], which_reads_aligned))
								correct_count += 1
							else:
								incorrect_count += 1 # GAGE_A
							#####
						#####

						# Reset
						curr_read_name = read_name
						curr_read_str = read_records[read_name]
						max_score = (0, "")
					# else, continue using this read
					#####
				else:
					# initialize
					curr_read_name = read_name
					curr_read_str = read_records[read_name]
				#####

				# Get the alignment region's kmers
				ref_k_set = getKmers(ref_record[ref_name][ref_start:ref_end], k_size)

				# score alignments with different weighting schemes
				shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum, shared_error_sum, non_shared_error_sum = counts(getKmers(curr_read_str[read_start:read_end], k_size), ref_k_set, kmer_table)
				#sys.stderr.write("%d, %d, %d, %d, %d, %d\n" % (shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum, shared_error_sum, non_shared_error_sum))
				score = weightJaccard(non_unique_kmer_weight, unique_kmer_weight, shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum, shared_error_sum, non_shared_error_sum)
				#sys.stderr.write("Max score = %0.8f to part %s, New Score = %0.8f to part %s\n" % (max_score[0], max_score[1], score, ref_name))
				if score > max_score[0]:
					max_score = (score, ref_name)
				#####
			#####


		
	#####


	# Calculate Performance
	
	# p_true = float(incorrect_count) / float(incorrect_count + correct_count) # turnover rate if which_reads_aligned == which error
	# p_false = float(correct_count) / float(incorrect_count + correct_count)		# retention rate if which_reads_aligned == which error

	#print("%0.8f\t%0.8f\t%0.8f\t%s\t%s\t%s\t" % (err_rate, p_true, p_false, which_error, which_reads_aligned, sam_file))

	# Or just return the counts, and calculate the rates later
	print("%0.8f\t%d\t%d\t%d\t%s\t%s\t%s\t" % (err_rate, correct_count, incorrect_count, unmapped_count, which_error, which_reads_aligned, sam_file))

	# Check later for if which_reads_aligned == which error in downstream analysis for True positives and true negatives





if __name__ == "__main__": main()
