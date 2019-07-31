# Script to create error

import sys
import random


alphabet = {0 : 'A',
			1 : 'T', 
			2 : 'C', 
			3 : 'G'}
num_alphabet = 4

err_types = { 0 : 'mismatch',
			  1 : 'insertion', 
			  2 : 'deletion'}
num_err_types = 3


def sim_err(error_rate, seq):

	new_seq = ""

	threshold = 100 * error_rate
	error_count = 0

	for c in seq:
		x = random.randint(0,101)

		if x <= threshold :
			# Do edit base
			error_count += 1

			# Pick a base
			x = random.randint(0,101)
			i = x % num_alphabet
			base = alphabet[i]

			# pick an error type
			x = random.randint(0,101)
			i = x % num_err_types
			err_type = err_types[i]


			if err_type == 'mismatch':
				new_seq += base
			elif err_type == "insertion":
				new_seq += base
				new_seq += c
			elif err_type == 'deletion':
				pass
			#####

		else:
			new_seq += c
		#####
	#####

	sys.stderr.write("Length: %d\nNumber of errors introduced: %d\n" % (len(seq), error_count))

	return new_seq


def main():

	# Get params
	fasta_file = sys.argv[1]
	error_rate = float(sys.argv[2])
	outfile = sys.argv[3]


		
	fh = open(outfile, "w")

	sys.stderr.write("Writing to %s\n" % outfile)
	
	curr_read = ""
	seq = ""
	for line in open(fasta_file, "r"):
		if ">" in line:
			if seq:
				# already started

				new_seq = sim_err(error_rate, seq)
				fh.write("%s\n%s\n" % (curr_read, new_seq))

			#####
			curr_read = line.strip() # header
			seq = ""
		else:
			seq += line.strip()
		#####
	#####

	# last one
	new_seq = sim_err(error_rate, seq)
	fh.write("%s\n%s\n" % (curr_read, new_seq))
	fh.close()




if __name__ == "__main__": main()