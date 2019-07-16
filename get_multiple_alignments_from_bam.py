# Script to extract reads with multiple alignments from a bam file and output an abbreviated 


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

	# Stdin for sam file
	sam_file = sys.stdin

	true_origin_bedfile = sys.argv[1] # bedfile for the true regions of the reads

	true_origin_table = {}
	for line in open(true_origin_bedfile, "r"):
		read = line.split()[-1]
		start = int(line.split()[1])
		end = int(line.split()[2])
		true_origin_table[read] = (start, end)
	#####

	curr_read_name = ""
	first_data = ""

	for line in sam_file:
		read_name = line.split()[0]

		start = int(line.split()[3])
		isCorrect = isTrue(start, true_origin_table[read_name][0], true_origin_table[read_name][1])

		data = "%s\t%s" % (line.strip(), isCorrect)

		if (read_name == curr_read_name):
			if (second):
				print(first_data)
				second = False
			#####
			print(data)
		else:
			curr_read_name = read_name
			first_data = data
			second = True
		#####
	#####

if __name__ == "__main__": main()