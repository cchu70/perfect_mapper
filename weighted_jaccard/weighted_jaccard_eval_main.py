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

	sch_test = {}
	curr_sch_scores = {}

	for line in open(wj_file, "r"):

		read_name, ground_truth, scores_table = parseWJ(line.strip())

		if not curr_read:
			curr_read = read_name
		elif curr_read != read_name:
			# get the current's results
			for scheme in curr_sch_scores:

				try:
					x = sch_test[scheme]
				except KeyError:
					sch_test[scheme] = [0,0,0,0]
				#####

				# All the scores based on current scheme from all the alignments for the current read
				align_scores = curr_sch_scores[scheme]

				max_align = max(align_scores)

				for a in align_scores:
					if a == max_align:
						if a[1] == "True":
							# increment TP
							sch_test[scheme][0] += 1
						else:
							# increment TN
							sch_test[scheme][1] += 1
							# print("TN")
							# print(curr_read)
							# print(curr_sch_scores[scheme])
						#####
					else:
						if a[1] == "True":
							# Increment FP
							sch_test[scheme][2] += 1
							# print("FP")
							# print(curr_read)
							# print(curr_sch_scores[scheme])
						else:
							# increment FN
							sch_test[scheme][3] += 1

						#####
					#####
				#####
			#####

			# update 
			curr_read = read_name
			curr_sch_scores = {}
		#####

		for scheme in scores_table:
			score = scores_table[scheme]
			try:
				curr_sch_scores[scheme].append((score, ground_truth))
			except:
				curr_sch_scores[scheme] = [(score, ground_truth)]
		#####
	#####


	# # Print results
	for sch in sch_test:
		rates = sch_test[sch]
		tp = rates[0]
		tn = rates[1]
		fp = rates[2]
		fn = rates[3]
		print("%s\t%d\t%d\t%d\t%d" % (sch, tp, tn, fp, fn))
	#####






if __name__ == "__main__": main()





