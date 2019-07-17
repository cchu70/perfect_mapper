# Script to evaulate the performance of the output from different weights used to calculate a jaccard similarity score


import sys

def parseWJ(wj_str):
	data = wj_str.split("\t")
	read_name = data[0]
	ground_truth = data[3]

	scores = data[4:]
	scores_table = {}
	for s in scores:
		scheme, score = s.split("=")
		scores_table[scheme] = score
	#####

	return read_name, ground_truth, scores_table

def main():


	wj_file = sys.argv[1]


	curr_read = None
	comp_table = {}
	sch_test_true = {}
	sch_test_false= {}
	for line in open(wj_file, "r"):

		read_name, ground_truth, scores_table = parseWJ(line.strip())

		if not curr_read:
			curr_read = read_name
		elif curr_read != read_name:
			# get the current's results
			for sch in comp_table:

				try:
					x = sch_test_true[sch]
				except KeyError:
					sch_test_true[sch] = 0
				#####

				try:
					x = sch_test_false[sch]
				except KeyError:
					sch_test_false[sch] = 0
				#####
		
				if comp_table[sch][1] == "True":
					sch_test_true[sch] += 1
				else:
					sch_test_false[sch] += 1
				#####
			#####

			# update 
			curr_read = read_name
			comp_table = {}
		#####


		# Get each score for each scheme
		for sch in scores_table:
			sch_score = (scores_table[sch], ground_truth)
			try:
				comp_table[sch] = max(comp_table[sch], sch_score)
			except KeyError:
				comp_table[sch] = sch_score
			#####
		#####
	#####


	# Print results
	for sch in sch_test_true:
		correct_count = float(sch_test_true[sch])
		incorrect_count = float(sch_test_false[sch])
		total = correct_count + incorrect_count

		tp_rate = correct_count / total
		fp_rate = incorrect_count / total

		print("%s\t%0.5f\t%0.5f" % (sch, tp_rate, fp_rate))
	#####






if __name__ == "__main__": main()





