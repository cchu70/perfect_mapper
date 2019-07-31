# script to take in all the out files from testing the turnover rate when introducing certain levels of error

import sys


def main():


	fofn = sys.argv[1] # file that lists the file paths of all the out files


	for file in open(fofn, "r"):

		for line in open(file.strip(), "r"):

			# Parse
			err, which_err, from_A_align_A, from_A_align_B, from_B_align_B, from_B_align_A, x, y = line.strip().split()

			data = line.strip().split()
			err = float(data[0])
			which_err = data[1] 
			from_A_align_A = int(data[2])
			from_A_align_B = int(data[3])
			from_B_align_B = int(data[4])
			from_B_align_A = int(data[5])



			if which_err == "A":
				p_turnover = float(from_A_align_B) / float(from_A_align_A + from_A_align_B)
				p_remaining = float(from_A_align_A) / float(from_A_align_A + from_A_align_B)
			else:
				p_turnover = float(from_B_align_A) / float(from_B_align_B + from_B_align_A)
				p_remaining = float(from_B_align_B) / float(from_B_align_B + from_B_align_A)
			#####

			print("%0.2f\t%0.2f\t%0.2f\t%s" % (err, p_turnover, p_remaining, which_err))


		#####
	#####










if __name__ == "__main__": main()