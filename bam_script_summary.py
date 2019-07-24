# script to evaluate a bam file

import sys

def main():

	# data to retrieve
	prim_count
	sec_count

	foward_align_count

	true_count

	false_count


	true_positive
	false_positive
	true_negative
	false_negative

	mult_align_count
	mult_align_count_read
	single_align_count


	# get the file from stdin

	sam = sys.stdin

	# bedfile compare

	bedfile = sys.argv[1]

	true_origin_table = {}
	for line in open(true_origin_bedfile, "r"):
		data = line.strip().split("\t")
		read = data[0]
		start = int(data[1])
		end = int(data[2])
		true_origin_table[read] = (start, end)
	#####


	for line in sam:

		# parse the flag