# Script to set a scoring scheme

from  weighted_jaccard_func import getKmers, parsePaf, Alignment, Scheme, parseUniqueFile, parseSam, parseFasta, counts, weightJaccard, align_file_parser
import sys
import os
import subprocess
import time
import threading



def parse_k_count_file(k_count_line):
	data = k_count_line.split("\t")
	read_name = data[0]
	map_truth = data[1]
	ref_start = int(data[2])
	ref_end = int(data[3])
	ground_truth = data[4]
	pid = float(data[5])
	shared_unique_sum = int(data[6])
	shared_non_unique_sum = int(data[7])
	non_shared_unique_sum = int(data[8])
	non_shared_non_unique_sum = int(data[9])
	return read_name, map_truth, ref_start, ref_end, ground_truth, pid, shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum

def main():


	k_count_file = sys.argv[1]
	sch_start = float(sys.argv[2])
	sch_end = float(sys.argv[3])
	step = float(sys.argv[4])
	# op = sys.argv[5]



	# Decimals
	schemes = [] 

	start = sch_start

	while start <= sch_end:
		j = sch_start
		while j <= sch_end:
			schemes.append(Scheme(j, start))
			j += step
		start += step
	#####



	curr_read_str = None
	curr_read_name = None
	# curr_read_k_set = {}

	for line in open(k_count_file, "r"):

		read_name, map_truth, ref_start, ref_end, ground_truth, pid, shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum = parse_k_count_file(line.strip())
		alignment = Alignment(read_name, map_truth, ref_start, ref_end, ground_truth, pid)

		# score alignments with different weighting schemes
		for sch in schemes:
			x = weightJaccard(sch.non_unique_weight, sch.unique_weight, shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum)

			if (x < 0):
				sys.stderr.write("%s\n" % alignment)
				sys.stderr.write("%d\t%d\t%d\t%d\n"% (shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum))
				assert False
			alignment.scores[sch] = x
			# alignment.scores[sch] = [shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum]
		#####
		print(alignment.scoreString())
	#####



if __name__ == "__main__": main()




