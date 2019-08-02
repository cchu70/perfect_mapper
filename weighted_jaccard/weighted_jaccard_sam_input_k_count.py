

#!/usr/bin/env python

# SCript to run through a bunch of files to weighted_jaccard_count_plain_sam_input.py

import threading
import time

import sys


def run_count(err_str, prefix):

	cmd = "python ../../scripts/mashmap_postfilter/weighted_jaccard/weighted_jaccard_count_plain_sam_input.py GAGE_A.sim_reads.fasta error_%s/%s.err_%s_A.v_1.fasta error_%s/%s_minimap2.N50_r3k.split.err_%s_A.v_1.aligned_A.sam GAGE.kmerlist.txt 21 GAGE_A GAGE_A %s" % (err_str, prefix, err_str, err_str, prefix, err_str, err_str)
	print(cmd)

def main():

	errors = ['0.0', '0.001', '0.002', '0.003', '0.004', '0.005', '0.006', '0.007', '0.008', '0.009', '0.01']
	prefix = 'AAF'
	v = 10


	for e in errors:
		i = 1
		while i <= v:	
			x = threading.Thread(target=run_count, args=(e, prefix))
			x.start()
			i += 1
		#####
	#####








if __name__ == "__main__": main()