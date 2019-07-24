# script to evaluate a bam file

import sys


def isTrue(align_start, true_start, true_end):
	length = true_end - true_start
	#print(align_start, true_start, true_end)
	# print ("%d < %d < %d ?" % ((true_start - length / 2), align_start, (true_start + length / 2)))
	if ( (true_start - length / 2) <= align_start < (true_start + length / 2) ):
		return True
	#####
	return False


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

	no_liftover
	unmapped = 0


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


	curr_read_name = ""

	for line in sam:
		data =  line.split()

		read_name =data[0]

		#sam file
		start = int(data[3])
		flag = int(data[1])

		try:
			# has liftover coordinates
			index = true_origin_table[read_name]
			is_correct = isTrue(start, true_origin_table[read_name][0], true_origin_table[read_name][1])
		except KeyError:
			if(curr_read_name != read_name)
				no_liftover += 1
			#####
			continue
		#####



		# parse the flag

		if flag == 0 or flag == 16:
			primary = True


		if flag == 4:
			unmapped += 1

		if not flag & 16:
			forward = False


		if flag & 2048:
			split = True


		if curr_read_name == read_name and on_second:
			multiple = True
		elif curr_read_name == read_name and not on_second:
			on_second = True
		else:
			single += 1 # the current read has only one alignment
		#####





		if(primary)


		if(is_correct)


		if(multiple)


		if(forward)


		if(split)







