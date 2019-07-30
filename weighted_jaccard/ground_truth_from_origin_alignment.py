# Script to assign ground truth based on the alignment of the non-simulated versions of the reads


import sys

def is_true(length, origin_start, align_start):
	if ( (origin_start - length / 2) <= align_start < (origin_start + length / 2) ):
		return True
	#####
	return False



def main():

	ground_truth_file = sys.argv[1]
	align_file = sys.stdin


	ground_truth = {}
	for line in open(ground_truth_file, "r"):
		read_name, start = line.strip().split()
		ground_truth[read_name] = int(start)
	#####


	for line in align_file:
		read_name = line.split()[0]

		chrm, read_start, read_end, fwd = read_name.split("_")
		length = int(read_end) - int(read_start)

		# sam file
		align_start = int(line.split()[3])
		flag = int(line.split()[1])

		try:
			index = ground_truth[read_name]
		except:
			continue 
		#####

		gt = is_true(length, ground_truth[read_name], int(align_start))

		print("%s\t%s" % (line.strip(), gt))

	#####



if __name__ == "__main__": main()