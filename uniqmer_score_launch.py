# Main pipeline for counting and scoring shared and non-shared uniqmers and non uniqmers
# called upon by uniqmer_score_main.py

import sys


def main():

	stage = sys.argv[1]
	config = sys.argv[2] 	# configuration file
	summary = sys.argv[3] 		# summary file 

	params = parse_config(config)


	if stage == 'uniqmer_db':
		# Make uniqmer database

		next_stage = 'alignment'

	elif stage == 'align_sim':
		# align the reads to the ref with chosen alignment (minimap2 or mashmap)



		next_stage = 'ground_truth'

	elif stage == 'ground_truth':
		# retrieve the ground truth for each read

		next_stage = 'mult_align_only'

	elif stage == "mult_align_only":
		# Or just a basic filter stage


		next_stage = 'count'

	elif stage == "count":
		# count the uniqmers and non-uniqmers, and shared and non-shared mers

	next_stage = 'scheme_score'

	elif stage == "scheme_score":
		# Use the kmer counts to score the weighted Jaccard similarity score

		next_stage = 'eval_scores'

	elif stage == "eval_scores":


		next_stage = 'calc_pid'

	elif stage == "calc_pid":
		# if minimap2, need to calculate the percent Identity

		next_stage = 'eval_pid'
	elif stage == "eval_pid"

	else:
		# not a stage
		sys.stderr.write("%s is not a vaild stage." % stage)
		assert False
	#####


	# print summary 


	# launch next stage

	os.system('python uniqmer_score_launch.py %s %s %s' % (options.stage, options.config_file, options.summary_file))








if __name__ == "__main__": main()