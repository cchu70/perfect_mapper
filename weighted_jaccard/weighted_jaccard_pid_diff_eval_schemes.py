# Script to get files for each scheme on the difference in pid between the best and second best alignment based on our weighted score


import sys


def main():

	# Get the file
	sch_score_file = sys.argv[1]


	curr_read = None

	scores = []


	for line in open(sch_score_file, "r"):
		read_name, ref_start, ref_end, ground_truth, pid, score = line.split() # awk out the columns


		score = float(score.split("=")[-1])

		alignment = Alignment(read_name, map_truth, ref_start, ref_end, ground_truth, pid)

		if not curr_read:
			curr_read = read_name
		else:
			if (read_name != curr_read):
				# get the max

				scores.sort()

				best = scores[0]
				second_best = scores[1]

				if best[1].ground_truth == "True":
						diff = best[0] - second_best[0]
						print("correct\t%0.8f" % diff)
				else:
					# sec must be true
					diff = second_best[0] - best[0]
					print("incorrect\t%0.8f" % diff)
				#####

				#reset
				curr_read = read_name
				scores = []
			#####


			scores.append((score, alignment))
		#####
	#####







if __name__ == "__main__": main()