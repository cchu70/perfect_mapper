

#!/usr/bin/env python

# SCript to run through a bunch of files to weighted_jaccard_count_plain_sam_input.py

import threading
import time
import subprocess

import sys


def run_count(err_str, prefix, ver, which_part, kmer_list):

	line = ""
	# TO_DO: remove dependency on GAGE and kmerlist file

	cmd = "python ../../scripts/mashmap_postfilter/weighted_jaccard/weighted_jaccard_count_plain_sam_input.py GAGE_%s.sim_reads.fasta error_%s/%s.err_%s_%s.v_%d.fasta error_%s/%s_minimap2.N50_r3k.split.err_%s_%s.v_%s.aligned_%s.sam %s 21 GAGE_%s GAGE_%s %s" % (which_part, err_str, prefix, err_str, which_part, ver, err_str, prefix, err_str, which_part, ver, which_part, kmer_list, which_part, which_part, err_str)
	#print(cmd)
	p1 = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

	output = p1.communicate()[0]
	print(output)
	return output


def main():


	errors = sys.argv[1]

	errors = errors.split(",")				# errors = ['0.0', '0.0001', '0.0002', '0.0003', '0.0004', '0.0005', '0.0006', '0.0007', '0.0008', '0.0009']

	v_start = int(sys.argv[2])				# There are multiple versions of a simulated error
	v_end = int(sys.argv[3])
	v = v_start	

	kmer_list = sys.argv[4]					# text file with the format "kmer_string"<tab>"frequency" where all frequencies must be >= 1

	prefix = sys.argv[5]					# prefix = 'AAG'							

	threads = list()

	for e in errors:
		i = 1
		while i <= v_end:	
			xA = threading.Thread(target=run_count, args=(e, prefix, i, 'A', kmer_list))
			threads.append(xA)

			xB = threading.Thread(target=run_count, args=(e, prefix, i, 'B', kmer_list))
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