

#!/usr/bin/env python

# SCript to run through a bunch of files to weighted_jaccard_count_plain_sam_input.py

import threading
import time
import subprocess

import sys


def run_count(err_str, prefix, ver):

	line = ""
	cmd = "python ../../scripts/mashmap_postfilter/weighted_jaccard/weighted_jaccard_count_plain_sam_input.py GAGE_A.sim_reads.fasta error_%s/%s.err_%s_A.v_%d.fasta error_%s/%s_minimap2.N50_r3k.split.err_%s_A.v_%s.aligned_A.sam GAGE.kmerlist.txt 21 GAGE_A GAGE_A %s" % (err_str, prefix, err_str, ver, err_str, prefix, err_str, ver, err_str)
	p1 = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

	output = p1.communicate()[0]
	print(output)
	return output


def main():

	errors = ['0.0', '0.0001', '0.0002', '0.0003', '0.0004', '0.0005', '0.0006', '0.0007', '0.0008', '0.0009', '0.001']
	prefix = 'AAF'
	v = 10

	threads = list()

	for e in errors:
		i = 1
		while i <= v:	
			x = threading.Thread(target=run_count, args=(e, prefix, i))
			threads.append(x)
			x.start()
			i += 1
		#####
	#####

	for index, thread in enumerate(threads):
		thread.join()
	#####








if __name__ == "__main__": main()