# script to parse mashmap to only get forward multiple alignments and get ground truth. Also mark which one mashmap would return as primary

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

	mashmap_out = sys.argv[1] # plain mashmap output
	true_origin_bedfile = sys.argv[2] # bed file of the original locations of the reads

	true_origin_table = {}
	for line in open(true_origin_bedfile, "r"):
		data = line.strip().split("\t")
		read = data[0]
		start = int(data[1])
		end = int(data[2])
		true_origin_table[read] = (start, end)
	#####


	curr_read = None
	curr_read_pid = []
	for line in open(mashmap_out, "r"):
		data = line.split()

		read_name = data[0]
		start = int(data[7])
		pid = float(data[-1])

		# print(line.strip())
		# print(start, true_origin_table[read_name][0], true_origin_table[read_name][1])
		ground_truth = isTrue(start, true_origin_table[read_name][0], true_origin_table[read_name][1])

		is_forward = False
		if data[4] == "+": 
			is_forward = True
		#####

		# only consider forward alignments
		if (data[4] == "+"):
			if curr_read != read_name:
				# Print out current read

				if len(curr_read_pid) > 1:
					max_align = max(curr_read_pid)

					for align in curr_read_pid:
						mash_best = "S"
						if align == max_align:
							mash_best = "P"
						#####

						print("%s\t%s\t%s" % (align[2], mash_best, align[1]))
					#####
				######

				# reset to new read
				curr_read = read_name
				curr_read_pid = []
			#####
			curr_read_pid.append((pid, ground_truth, line.strip()))
		#####










if __name__ == "__main__": main()