

#!/usr/bin/env python

# SCript to run through a bunch of files to weighted_jaccard_count_plain_sam_input.py

import threading
import time
import subprocess

import sys


def run_count(err_str, prefix, ver, which_part_aligned, which_part_errored, uniq_weight, non_uniq_weight):

	script = "python ../../scripts/mashmap_postfilter/weighted_jaccard/weighted_jaccard_count_plain_sam_input.py"
	sim_reads = "GAGE_%s.sim_reads.fasta" % (which_part_aligned)
	target = "error_%s/%s_split.err_%s_%s.v_%d.fasta" % (err_str, prefix, err_str, which_part_errored, ver)
	align_file = "error_%s/%s_minimap2.N50_r3k.split.err_%s_%s.v_%s.aligned_%s.sam" % (err_str, prefix, err_str, which_part_errored, ver, which_part_aligned)
	kmerlist = "GAGE.kmerlist.txt"
	kmer_size = 21
	which_part_aligned = "GAGE_%s" % (which_part_aligned)
	which_part_errored = "GAGE_%s" % (which_part_errored)

	cmd = "%s %s %s %s %s %s %s %s %s %s %s" % (script, sim_reads, target, align_file, kmerlist, kmer_size, which_part_aligned, which_part_errored, err_str, uniq_weight, non_uniq_weight)
	# cmd = "python ../../scripts/mashmap_postfilter/weighted_jaccard/weighted_jaccard_count_plain_sam_input.py GAGE_%s.sim_reads.fasta error_%s/%s_split.err_%s_%s.v_%d.fasta error_%s/%s_minimap2.N50_r3k.split.err_%s_%s.v_%s.aligned_%s.sam GAGE.kmerlist.txt 21 GAGE_%s GAGE_%s %s" % (which_part, err_str, prefix, err_str, which_part, ver, err_str, prefix, err_str, which_part, ver, which_part, which_part, which_part, err_str)
	# print(cmd)
	p1 = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

	output = p1.communicate()[0]
	print(output)
	# return output


def main():


	errors = sys.argv[1]

	errors = errors.split(",")

	prefix = sys.argv[2]

	uniq_weights = sys.argv[3]
	uniq_weights = uniq_weights.split(",")

	non_uniq_weights = sys.argv[4]
	non_uniq_weights = non_uniq_weights.split(",")


	schemes = []

	for uniq_w in uniq_weights:
		for non_uniq_w in non_uniq_weights:
			schemes.append("%s:%s" % (uniq_w, non_uniq_w))
		#####
	####

	schemes_arg = ",".join(schemes)
	print(schemes_arg)
	assert False
	# errors = ['0.0', '0.0001', '0.0002', '0.0003', '0.0004', '0.0005', '0.0006', '0.0007', '0.0008', '0.0009']
	# errors = ['0.0', '0.0001']
	# prefix = 'AAG'
	v = 10

	threads = list()

	parts = ['A', 'B']

	for e in errors:
		i = 1
		while i <= v:

			for part in parts:	
				xA = threading.Thread(target=run_count, args=(e, prefix, i, 'A', part))
				threads.append(xA)

				xB = threading.Thread(target=run_count, args=(e, prefix, i, 'B', part))
				threads.append(xB)
			#####
			i += 1
		#####
	#####


	for index, thread in enumerate(threads):
		thread.start()
	#####

	for index, thread in enumerate(threads):
		thread.join()
	#####








if __name__ == "__main__": main()