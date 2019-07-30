# Script to create error

import sys
import random


def sim_err(error_rate, seq):

	new_seq = ""

	threshold = 100 * error_rate


	error_count = 0

	for c in seq:
		x = random.randint(0,101)

		if x <= threshold :
			# Do edit base

			# Pick a base

			# pick an error type
		else:
			new_seq += c


def main():

	# Get params
	fasta_file = sys.argv[1]
	error_rate = float(sys.argv[2])
	iterations = int(sys.argv[3])
	prefix = sys.argv[4]


	for i in range(iteratios):
		new_fasta_file = "%s.%0.2f.v_%d.fa" % (prefix, error_rate, i)
		fh = open(new_fasta_filem, "w")
		curr_read = ""
		seq = ""
		for line in open(fasta_file, "r"):
			if ">" in line:
				if seq:
					# already started

					new_seq = sim_err(error_rate, seq)
					fh.write("%s\n%s\n" % (curr_read, seq))

				#####
				curr_read = line.strip() # header
				seq = ""
			else:
				seq += line.strip()
			#####
		#####

		# last one
		new_seq = sim_err(error_rate, seq)
		fh.write("%s\n%s\n" % (curr_read, seq))
		fh.close()
	#####



if __name__ == "__main__": main()