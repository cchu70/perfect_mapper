# Main script to run the analysis of using weighted jaccard scoring on uniq-mers

import sys
import os
from weighted_jaccard_lib import init_parser
import subprocess

def main():

	# Parse arguments for stage
	parser = init_parser()
	(options, args) = parser.parse_args()

	# parse the configuration file

	params = parse_config(options.config_file)

	# run main pipeline
	os.system('python uniqmer_score_launch.py %s %s %s' % (options.stage, options.config_file, options.summary_file))

#####





if __name__ == "__main__": main()