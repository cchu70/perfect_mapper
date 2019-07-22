# Script to evaulate the performance of the output from different weights used to calculate a jaccard similarity score

# python this.script.py weighted_jaccard.scheme_scores.txt

import sys
from  weighted_jaccard_func import Alignment

def parseWJ(wj_str):
	data = wj_str.split("\t")
	read_name = data[0]
	map_truth = data[1]
	start = int(data[2])
	end = int(data[3])
	ground_truth = data[4]

	scores = data[5:-1] # do not include the weights
	scores_table = {}
	for s in scores:
		# scheme, score, equation = s.split("=")
		# scores_table[scheme] = (float(score), equation)
		scores_table[s.split("=")[0]] = s.split("=")[1]

	return read_name, map_truth, start, end, ground_truth, scores_table

def main():


	wj_file = sys.argv[1]


	curr_read = None

	sch_test = {}
	curr_sch_scores = {}

	for line in open(wj_file, "r"):

		read_name, map_truth, ref_start, ref_end, ground_truth, scores_table = parseWJ(line.strip())
		alignment = Alignment(read_name, map_truth, ref_start, ref_end, ground_truth)

		if not curr_read:
			curr_read = read_name
		elif curr_read != read_name:
			# get the current's results
			wrong = False
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
						if a[1].ground_truth == "True":
							# increment TP
							sch_test[scheme][0] += 1
						else:
							# increment FP
							sch_test[scheme][1] += 1
							#print("%s\t%s\t%s\t%s\tTN" % (curr_read, scheme, a[0], a[1]))
							# print(a[0])
							# print(a[1])
							wrong = True
						#####
					else:
						if a[1].ground_truth == "True":
							# Increment FN
							sch_test[scheme][2] += 1
							#print("%s\t%s\t%s\t%s\tFP" % (curr_read, scheme, a[0], a[1]))
							# print(curr_read)
							# print(curr_sch_scores[scheme])
							wrong = True
						else:
							# increment TN
							sch_test[scheme][3] += 1

						#####
					#####
				#####
			#####
			if(wrong):
				print("%s\n" % curr_read)

			# update 
			curr_read = read_name
			curr_sch_scores = {}
		#####

		for scheme in scores_table:
			score = scores_table[scheme]
			try:
				curr_sch_scores[scheme].append((score, alignment))
			except:
				curr_sch_scores[scheme] = [(score, alignment)]
		#####
	#####


	#assert False
	# # Print results
	for sch in sch_test:
		rates = sch_test[sch]
		tp = rates[0]
		fp = rates[1]
		fn = rates[2]
		tn = rates[3]
		sys.stderr.write("%s\t%d\t%d\t%d\t%d" % (sch, tp, fp, fn, tn))
	#####

# chr22_part05_37850_40445
# chr22_part05_10344_11445
# chr22_part05_34922_36159
# chr22_part22_182338_187995
# chr22_part24_20962_22189
# chr22_part25_2787121_2789117




if __name__ == "__main__": main()





