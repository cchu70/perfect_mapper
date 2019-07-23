# Script to get files for each scheme on the difference in pid between the best and second best alignment based on our weighted score


import sys
from weighted_jaccard_func import parseSam, Alignment


def main():

	# Get the file
	sch_score_file = sys.argv[1]


	curr_read = None

	curr_best = None
	curr_best_score = -1


	curr_second_best = None

	curr_second_best_score = -1
	
	curr_true = None

	scores = []


	for line in open(sch_score_file, "r"):
		read_name, map_truth, ref_start, ref_end, ground_truth, pid, score = line.split("\t") # awk out the columns


		score = float(score.split("=")[-1])

		# alignment = Alignment(read_name, "S", ref_start, ref_end, ground_truth, pid)

		if not curr_read:
			curr_read = read_name
		else:
			if (read_name != curr_read):
				# get the max
				scores.sort()

				best_align = scores[0]
				second_best_align = scores[1]

				if best_align[5] == "True":
					diff = best_align[6] - second_best_align[6]
					print("correct\t%0.8f" % diff)
				else:
					diff = second_best_align[6] - best_align[6]
					print("incorrect\t%0.8f" % diff)
				#####

				#reset
				curr_read = read_name
				scores = []
			#####
			
			
			scores.append((score, read_name, map_truth, ref_start, ref_end, ground_truth, pid))
		#####
	#####







if __name__ == "__main__": main()