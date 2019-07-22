# Script to set a scoring scheme

from  weighted_jaccard_func import getKmers, parsePaf, Alignment, Scheme, parseUniqueFile, parseSam, parseFasta, counts, weightJaccard, align_file_parser
import sys
import os
import subprocess
import time
import threading
import numpy



def parse_k_count_file(k_count_line):
	data = k_count_line.split("\t")
	read_name = data[0]
	map_truth = data[1]
	ref_start = int(data[2])
	ref_end = int(data[3])
	ground_truth = data[4]
	shared_unique_sum = data[5]
	shared_non_unique_sum = data[6]
	non_shared_unique_sum = data[7]
	non_shared_non_unique_sum = data[8]
	return read_name, map_truth, ref_start, ref_end, ground_truth, shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum

def main():


	k_count_file = sys.argv[1]
	sch_start = int(sys.argv[2])
	sch_end = int(sys.argv[3])
	step = float(sys.argv[4])
	op = sys.argv[5]



	# Decimals
	schemes = [Scheme(1, i) for i in numpy.arange(sch_start, sch_end, step)] # Linear
	# schemes = [Scheme(1, i) for i in range(sch_start, sch_end)] # Linear


	curr_read_str = None
	curr_read_name = None
	# curr_read_k_set = {}

	for line in open(k_count_file, "r"):

		read_name, map_truth, ref_start, ref_end, ground_truth, shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum = parse_k_count_file(line.strip())
		alignment = Alignment(read_name, map_truth, ref_start, ref_end, ground_truth)

		# score alignments with different weighting schemes
		for sch in schemes:
			x = weightJaccard(sch.non_unique_weight, sch.unique_weight, shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum)
			alignment.scores[sch] = x
			# alignment.scores[sch] = [shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum]
		#####
		print("%s\t%s" % (alignment.toString(), [shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum]))
	#####



if __name__ == "__main__": main()




