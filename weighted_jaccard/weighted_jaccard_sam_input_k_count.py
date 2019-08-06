

#!/usr/bin/env python

# SCript to run through a bunch of files to weighted_jaccard_count_plain_sam_input.py

import threading
import time
import subprocess

import sys


def run_count(err_str, prefix, ver, which_part):

	line = ""
	cmd = "python ../../scripts/mashmap_postfilter/weighted_jaccard/weighted_jaccard_count_plain_sam_input.py GAGE_%s.sim_reads.fasta error_%s/%s_split.err_%s_%s.v_%d.fasta error_%s/%s_minimap2.N50_r3k.split.err_%s_%s.v_%s.aligned_%s.sam GAGE.kmerlist.txt 21 GAGE_%s GAGE_%s %s" % (which_part, err_str, prefix, err_str, which_part, ver, err_str, prefix, err_str, which_part, ver, which_part, which_part, which_part, err_str)
	p1 = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

	output = p1.communicate()[0]
	print("%s\n%s" % (output,cmd))
	return output


def main():


	errors = sys.argv[1]

	errors = errors.split(",")

	prefix = sys.argv[2]

	# errors = ['0.0', '0.0001', '0.0002', '0.0003', '0.0004', '0.0005', '0.0006', '0.0007', '0.0008', '0.0009']
	# errors = ['0.0', '0.0001']
	# prefix = 'AAG'
	v = 10

	threads = list()

	for e in errors:
		i = 1
		while i <= v:	
			xA = threading.Thread(target=run_count, args=(e, prefix, i, 'A'))
			threads.append(xA)

			xB = threading.Thread(target=run_count, args=(e, prefix, i, 'B'))
			threads.append(xB)

			xA.start()
			xB.start()
			i += 1
		#####
	#####

	for index, thread in enumerate(threads):
		thread.join()
	#####








if __name__ == "__main__": main()