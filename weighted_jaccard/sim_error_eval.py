# script to take in all the out files from testing the turnover rate when introducing certain levels of error

import sys



def main():


	fofn = sys.argv[1] # file that lists the file paths of all the out files


	for file in open(fofn, "r"):
		sys.stderr.write("%s" % file)

		for line in open(file.strip(), "r"):

			# Parse
			# err, which_err, from_A_align_A, from_A_align_B, from_B_align_B, from_B_align_A, from_A_unaligned, from_B_unaligned, x, y = line.strip().split()

			data = line.strip().split()
			err = float(data[0])
			which_err = data[1] 
			from_A_align_A = int(data[2])
			from_A_align_B = int(data[3])
			from_B_align_B = int(data[4])
			from_B_align_A = int(data[5])
			# from_A_unaligned = int(data[6])
			# from_B_unaligned = int(data[7])

			A_correct = float(from_A_align_A)/ float(from_A_align_A + from_A_align_B)
			B_correct = float(from_B_align_B)/float(from_B_align_B + from_B_align_A)


			to_print_A = "%0.8f\t%0.8f\t%0.8f\tGAGE_%s\t%s" % (err, A_correct, 0, which_err, "GAGE_A")
			to_print_B = "%0.8f\t%0.8f\t%0.8f\tGAGE_%s\t%s" % (err, B_correct, 0, which_err, "GAGE_B")
			print(to_print_A)
			print(to_print_B)


			# if which_err == "A":
			# 	p_turnover = float(from_A_align_B) / float(from_A_align_A + from_A_align_B)
			# 	p_remaining = float(from_A_align_A) / float(from_A_align_A + from_A_align_B)
			# 	# p_total_wrong = (float(from_A_align_B) + float(from_A_unaligned)) / float(from_A_align_A + from_A_align_B)
			# else:
			# 	p_turnover = float(from_B_align_A) / float(from_B_align_B + from_B_align_A)
			# 	p_remaining = float(from_B_align_B) / float(from_B_align_B + from_B_align_A)
			# 	# p_total_wrong = (float(from_B_align_B) + float(from_B_unaligned)) / float(from_B_align_B + from_B_align_A)
			# #####


			# print("%0.8f\t%0.8f\t%0.8f\t%0.8f\t%s" % (err, p_turnover, p_remaining, p_total_wrong, which_err))
			# print("%0.8f\t%0.8f\t%0.8f\t%s" % (err, p_turnover, p_remaining, which_err))


		#####
	#####










if __name__ == "__main__": main()