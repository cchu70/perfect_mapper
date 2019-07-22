# Script to extract reads with multiple alignments from a bam file and output an abbreviated 

# samtools view bam_file | python this.script.py true_read.coordinates.bed
	# read ID in the 4th column

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
	paf_file = sys.stdin

	true_origin_bedfile = sys.argv[1] # bedfile for the true regions of the reads

	true_origin_table = {}
	for line in open(true_origin_bedfile, "r"):
		data = line.split()
		read = data[-1]
		start = int(data[1])
		end = int(data[2])
		true_origin_table[read] = (start, end)
	#####

	curr_read_name = ""
	first_data = ""

	for line in paf_file:
		read_name = line.split()[0]

		# paf file
		# start = int(line.split()[7])

		#sam file
		start = int(line.split()[3])
		flag = int(line.split()[1])

		try:
			index = true_origin_table[read_name]
		except:
			continue 
		#####

		isCorrect = isTrue(start, true_origin_table[read_name][0], true_origin_table[read_name][1])

		# ASSUMPTION THAT ALL THE ORIGINAL SIMULATED READS ARE FORWARD
		if flag & 16:
			isCorrect = False
		#####

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


# meryl-lookup -dump -threads 16 -memory 24G -sequence chrX.unpolished.fasta -mers chrX.unpolished.unique.meryl | awk '{if($3 == "T"){if($5 > 0){print $4"\t"$2}}}' > chrX.unpolished.asm_dump.txt ^C
# meryl count k=21 chrX.unpolished.fasta output chrX.unpolished.meryl -memory 24G -threads 16
# meryl equal-to 1 chrX.unpolished.meryl output chrX.unpolished.unique.meryl