# Script to get the percent ID of each alignment in a BAM file using the cigar string
# The compare the difference in the pid between the best and second best and whether it was right or wrong

# correct/incorrect pid_diff

import sys
from weighted_jaccard_func import parseSam, Alignment

def calc_align_length(cigar_str):
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
			#####
			num_string = ""
		#####
	#####

	return align_length


def main():

	sam_file = sys.stdin # pipe in sam file

	for line in sam_file:
		cigar_str = line.split()[5]
		num_mis_matches = line.split()[11]

		num_mis_matches = float(line.split()[11].split(":")[-1])
		align_length = calc_align_length(cigar_str)

		pid = 1 - float(num_mis_matches)/float(align_length)

		read_name, length, ref_start, ref_end, ground_truth, read_start, read_end, map_truth = parseSam(line.strip())

		alignment = Alignment(read_name, map_truth, ref_start, ref_end, ground_truth)

		print("%s\t%0.5f" % (alignment, pid))
		assert False





if __name__ == "__main__": main()