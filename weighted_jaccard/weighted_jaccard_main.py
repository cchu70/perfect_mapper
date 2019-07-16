# This is the main script to test different weighting schemes to compute a weighted jaccard score on shared kmers between a read and it's alignment region

from  weighted_jaccard_func import getKmers, score, parsePaf, Read, Alignment
from Bio import SeqIO
import sys
import os
import subprocess

def main():

	# Take inputs
	# 1) Alignment file
	# 2) Unique kmers 

	# Set weight schemes?
	sch_start = 1
	sch_end = 5
	k_size = 21
	read_fasta = "/data/Phillippy/projects/perfect-polish/chr22_info/chr22.sim_reads.fasta"
	ref_fasta = "/data/Phillippy/projects/perfect-polish/chr22_info/chr22.fasta"
	align_file = "/data/Phillippy/projects/perfect-polish/chr22_info/chr22.minimap2_N50_30kb.cigar.multiple_align_only.paf"

	if not os.path.isfile(read_fasta):
		assert False 

	if not os.path.isfile(ref_fasta):
		assert False

	#####


	schemes = [(1, w2) for w2 in range(sch_start, sch_end)]

	
	# Get read sequences
	sys.stderr.write("Parsing Read fasta: %s\n" % read_fasta)
	read_records = SeqIO.to_dict(SeqIO.parse(read_fasta, "fasta"))

	sys.stderr.write("Parsing Ref fasta: %s\n" % ref_fasta)
	ref_record = list(SeqIO.parse(ref_fasta, "fasta"))[0] # Should only be one reference

	curr_read = None

	for line in open(align_file, "r"):

		read_name, start, end, ground_truth = parsePaf(line.strip())

		sys.stderr.write("read name: %s, start: %d, end: %d, truth: %s\n" % (read_name, start, end, ground_truth))

		alignment = Alignment(start, end, ground_truth)

		# Check which read (current or next) this alignment corresponds to 

		if (curr_read):
			if (read_name != curr_read.read_name):
				# evaluate the curr read performance
				curr_read.print_alignments()
				assert False
			#####
		#####
		# Update to next read
		print(read_records[read_name].sequence)
		assert False
		curr_read = Read(read_name, getKmers(read_records[read_name].sequence, k_size))

		# Continue adding more alignments

		# Get the alignment region's kmers
		ref_k_set = getKmers(ref_record.seq[start:end], k_size)

		# score alignments with different weighting schemes
		for sch in schemes:
			alignment.add_score(sch, score(curr_read.k_set, ref_k_set, sch, k_size))
		#####

		curr_read.add_alignment(alignment)


		# Once I've seen and scored all of the current read's alignments, pick the best one, evaluate performace, etc.
	#####


	# For each read, pick the one with the highest score to be primary for each weighting scheme


	# Evaluate true and false rates












if __name__ == "__main__": main()




