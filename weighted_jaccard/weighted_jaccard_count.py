# Script to count the shared and unique kmers between two sequences# This is the main script to test different weighting schemes to compute a weighted jaccard score on shared kmers between a read and it's alignment region

from  weighted_jaccard_func import getKmers, parsePaf, Alignment, Scheme, parseUniqueFile, parseSam, parseFasta, counts, weightJaccard, align_file_parser
import sys
import os
import subprocess
import time
import threading
import numpy



def main():


	read_fasta = "/data/Phillippy/projects/perfect-polish/chr22_info/chr22.sim_reads.fasta"
	ref_fasta = "/data/Phillippy/projects/perfect-polish/chr22_info/chr22.fasta"
	align_file = "/data/Phillippy/projects/perfect-polish/chr22_info/representative_only.multiple_aligns_only.rev_false.minimap2_N50_30kb.real.sam"
	unique_k_file = "/data/Phillippy/projects/perfect-polish/chr22_info/chr22.asm.sck_list.txt"

	align_file_type = "sam"
	k_size = 21

	try:
		read_fasta = sys.argv[1]
		ref_fasta = sys.argv[2]
		align_file = sys.argv[3]
		align_file_type = sys.argv[4]
		unique_k_file = sys.argv[5]
		k_size = int(sys.argv[6])
	except:
		pass
	#####

	# Get read sequences
	sys.stderr.write("Parsing Read fasta: %s\n" % read_fasta)
	read_records = parseFasta(read_fasta) # Dictionary of read names and it's corresponding sequence

	sys.stderr.write("Parsing Ref fasta: %s\n" % ref_fasta)
	ref_record = list(parseFasta(ref_fasta).values())[0] # Should only be one reference

	sys.stderr.write("Parsing unique file: %s\n" % unique_k_file)
	unique_table = parseUniqueFile(unique_k_file)
	sys.stderr.write("Number of unique kmers: %d\n" % len(unique_table))

	curr_read_str = None
	curr_read_name = None
	# curr_read_k_set = {}

	parse_align_file = align_file_parser[align_file_type]

	for line in open(align_file, "r"):

		read_name, length, ref_start, ref_end, ground_truth, read_start, read_end, map_truth = parse_align_file(line.strip())

		# Gettin the pid
		num_mis_matches = float(line.split()[11].split(":")[-1])
		pid = 1 - float(num_mis_matches)/float(length)

		alignment = Alignment(read_name, map_truth, ref_start, ref_end, ground_truth, pid)


		# Check which read (current or next) this alignment corresponds to 
		if (curr_read_str):
			if (read_name != curr_read_name):
				# evaluate the curr read performance
				curr_read_name = read_name
				curr_read_str = read_records[read_name]
			# else, continue using this read
			#####
		else:
			# initialize
			curr_read_name = read_name
			curr_read_str = read_records[read_name]
		#####

		# Get the alignment region's kmers
		ref_k_set = getKmers(ref_record[ref_start:ref_end], k_size)

		# score alignments with different weighting schemes
		shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum = counts(getKmers(curr_read_str[read_start:read_end], k_size), ref_k_set, unique_table)
		print("%s\t%d\t%d\t%d\t%d" % (alignment, shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum))
		
	#####



if __name__ == "__main__": main()
