# Script to extract reads with multiple alignments from a bam file and output an abbreviated 


import sys

		
def parseSam(sam_string):
	# Based on sam file
	read_name = sam_string[0]
	flag = int(sam_string[1])
	start = int(sam_string[3])
	MQ = int(sam_string[4])
	cigar = sam_string[5]

	return read_name, flag, start, MQ, cigar

def main():

	# Stdin for sam file
	sam_file = sys.stdin

	curr_read_name = ""
	first_data = ""

	for line in sam_file:
		read_name = line.split()[0]
		if (read_name == curr_read_name):
			if (second):
				print(first_data)
				second = False
			#####
			print(line.strip())
		else:
			curr_read_name = read_name
			first_data = line.strip()
			second = True
		#####
	#####

if __name__ == "__main__": main()