# Script to filter out alignments in which a read is only assigned to a single region in the reference. This is used to evaluate the performance of the weighted Jaccard to be able to pick out the correct alignment, which would be uninformative if there is only one alignment to pick from to begin with
# prints to std.out

import sys


def main():

	# Stdin for sam file
	map_file = sys.argv[1]

	curr_read_name = ""
	first_data = ""

	for line in map_file:
		read_name = line.split()[0]
		print(read_name)

		if (read_name == curr_read_name):
			if (second):
				print(first_data)
				second = False
			#####
			print(data)
		else:
			curr_read_name = read_name
			first_data = line.strip()
			second = True
		#####
	#####




if __name__ == "__main__": main()