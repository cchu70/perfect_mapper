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

	ref_end = ref_start + length

	return read_name, length, ref_name, ref_start, ref_end, read_start, read_end


k_size = 21


def main():

	read_fasta = sys.argv[1]
	ref_fasta = sys.argv[2]
	sam_file = sys.argv[3]
	k_file = sys.argv[4]
	k_size = int(sys.argv[5])

	which_reads_aligned = sys.argv[6]
	which_error = sys.argv[7]
	err_rate = float(sys.argv[8])
	# only care if these are the same

	if which_reads_aligned != which_error:
		# Only consider the events where the reads we are mapping originates from the region we introduced error
		return
	#####

	unique_kmer_weight = 10 ###########################
	non_unique_kmer_weight = 1


	# Get read sequences
	sys.stderr.write("Parsing Read fasta: %s\n" % read_fasta)
	read_records = parseFasta(read_fasta) # Dictionary of read names and it's corresponding sequence

	sys.stderr.write("Parsing Ref fasta: %s\n" % ref_fasta)
	ref_record = list(parseFasta(ref_fasta).values())[0] # Should only be one reference

	sys.stderr.write("Parsing kmer file: %s\n" % k_file)
	kmer_table = parseKmerFile(k_file)
	sys.stderr.write("Finished counts kmers")


	# Initialize Loop
	curr_read_str = None
	curr_read_name = None


	max_score = (0, "")

	correct_count = 0
	incorrect_count = 0

	for line in open(sam_file, "r"):

		if "@" not in line:

			read_name, length, ref_name, ref_start, ref_end, read_start, read_end =  parseSam(line.strip())

			# Check which read (current or next) this alignment corresponds to 
			if (curr_read_str):
				if (read_name != curr_read_name):
					# evaluate the curr read performance
					# print("%s" % max_score[1])

					if max_score[1] == which_error: # GAGE_B
						correct_count += 1
					else:
						incorrect_count += 1 # GAGE_A

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
			shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum, shared_error_sum, non_shared_error_sum = counts(getKmers(curr_read_str[read_start:read_end], k_size), ref_k_set, kmer_table)

			score = weightJaccard(non_unique_kmer_weight, unique_kmer_weight, shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum)
			if score > max_score[0]:
				max_score = (score, ref_name)


		
	#####


	# end
	# print("%s" % max_score[1])
	p_turnover = float(incorrect_count) / float(incorrect_count + correct_count)
	p_remaining = float(correct_count) / float(incorrect_count + correct_count)

	print("%0.2f\t%0.2f\t%0.2f\t%s" % (err_rate, p_turnover, p_remaining, which_err))





if __name__ == "__main__": main()
