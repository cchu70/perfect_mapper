# Script to get the percent ID of each alignment in a BAM file using the cigar string
# The compare the difference in the pid between the best and second best and whether it was right or wrong

# correct/incorrect pid_diff

import sys
from weighted_jaccard_func import parseSam, Alignment

def calcPID(cigar_str):
	num_matches = 0
	align_length = 0
	num_string = ""

	for c in cigar_str:
		if c.isdigit():
			num_string += c
		else:
			# A letter
			d = int(num_string) # Convert the current num_string into a number

			if (c == "M" or c == "D" or c == "N"):
				# Only add up matches and deletions in the read
				align_length += d
				if (c == "M"):
					num_matches += 1
				#####
			#####
			num_string = ""
		#####
	#####
	return float(num_matches)/float(align_length)


def main():

	sam_file = sys.stdin # pipe in sam file

	for line in sam_file:
		cigar_str = line.split()[5]
		pid = calcPID(cigar_str)

		read_name, length, ref_start, ref_end, ground_truth, read_start, read_end, map_truth = parseSam(line.strip())

		alignment = Alignment(read_name, map_truth, ref_start, ref_end, ground_truth)

		print("%s\t%0.5f" % (alignment, pid))





if __name__ == "__main__": main()