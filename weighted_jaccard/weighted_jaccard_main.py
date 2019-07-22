# This is the main script to test different weighting schemes to compute a weighted jaccard score on shared kmers between a read and it's alignment region
# python this.script read_fasta ref_fasta align_file align_file_type list_of_unique_kmers scheme_start scheme_end k_size

from  weighted_jaccard_func import getKmers, parsePaf, Alignment, Scheme, parseUniqueFile, parseSam, parseFasta, counts, weightJaccard, align_file_parser
import sys
import os
import subprocess
import time
import threading



def main():


	read_fasta = "/data/Phillippy/projects/perfect-polish/chr22_info/chr22.sim_reads.fasta"
	ref_fasta = "/data/Phillippy/projects/perfect-polish/chr22_info/chr22.fasta"
	align_file = "/data/Phillippy/projects/perfect-polish/chr22_info/representative_only.multiple_aligns_only.rev_false.minimap2_N50_30kb.real.sam"
	unique_k_file = "/data/Phillippy/projects/perfect-polish/chr22_info/chr22.asm.sck_list.txt"

	align_file_type = "sam"

	sch_start = 1
	sch_end = 10
	k_size = 21

	try:
		read_fasta = sys.argv[1]
		ref_fasta = sys.argv[2]
		align_file = sys.argv[3]
		align_file_type = sys.argv[4]
		unique_k_file = sys.argv[5]
		sch_start = int(sys.argv[6])
		sch_end = int(sys.argv[7])
		k_size = int(sys.argv[8])
	except:
		pass
	#####

	# Set up schemes
	w = 2
	schemes = [Scheme(1, w ** i) for i in range(sch_start, sch_end)] # power

	# Decimals
	# schemes = [Scheme(1, w + i/10.0) for i in range(sch_start, sch_end)] # Linear

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

		read_name, length, ref_start, ref_end, ground_truth, read_start, read_end = parse_align_file(line.strip())

		alignment = Alignment(ref_start, ref_end, ground_truth, line.strip())

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
		for sch in schemes:
			#x = weightJaccard(sch.non_unique_weight, sch.unique_weight, shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum)
			# alignment.scores[sch] = x
			alignment.scores[sch] = [shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum]
		#####
		print(alignment.toString())
	#####



if __name__ == "__main__": main()




