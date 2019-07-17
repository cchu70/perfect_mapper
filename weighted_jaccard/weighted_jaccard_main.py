# This is the main script to test different weighting schemes to compute a weighted jaccard score on shared kmers between a read and it's alignment region

from  weighted_jaccard_func import getKmers, score, parsePaf, Read, Alignment, unique_table, parseUniqueFile
#from Bio import SeqIO
import sys
import os
import subprocess

def main():

	# Take inputs
	# 1) Alignment file
	# 2) Unique kmers 

	# Set weight schemes?
	sch_start = 0
	sch_end = 16
	k_size = 21
	read_fasta = "/data/Phillippy/projects/perfect-polish/chr22_info/chr22.sim_reads.fasta"
	ref_fasta = "/data/Phillippy/projects/perfect-polish/chr22_info/chr22.fasta"
	align_file = "/data/Phillippy/projects/perfect-polish/chr22_info/chr22.minimap2_N50_30kb.cigar.multiple_align_only.paf"
	unique_k_file = "/data/Phillippy/projects/perfect-polish/chr22_info/chr22.asm.sck_list.txt"


	w = 2
	schemes = [(1, w ** i) for i in range(sch_start, sch_end)]

	sys.stderr.write("Parsing unique file: %s\n" % read_fasta)
	unique_table = parseUniqueFile(unique_k_file)

	
	# Get read sequences
	sys.stderr.write("Parsing Read fasta: %s\n" % read_fasta)
	read_records = SeqIO.to_dict(SeqIO.parse(read_fasta, "fasta"))

	sys.stderr.write("Parsing Ref fasta: %s\n" % ref_fasta)
	ref_record = list(SeqIO.parse(ref_fasta, "fasta"))[0] # Should only be one reference

	curr_read = None
	alignments = []

	for line in open(align_file, "r"):

		read_name, length, start, end, ground_truth = parsePaf(line.strip())

		# sys.stderr.write("read name: %s, start: %d, end: %d, truth: %s\n" % (read_name, start, end, ground_truth))

		alignment = Alignment(start, end, ground_truth)


		# Check which read (current or next) this alignment corresponds to 
		if (curr_read):
			if (read_name != curr_read.read_name):
				# evaluate the curr read performance
				# curr_read.print_alignments()
				# assert False
				curr_read = Read(read_name, length, getKmers(read_records[read_name].seq, k_size))
			#####
		else:
			# initialize
			# sys.stderr.write("Initialize\nread name: %s, start: %d, end: %d, truth: %s\n" % (read_name, start, end, ground_truth))
			curr_read = Read(read_name, length, getKmers(read_records[read_name].seq, k_size))
		#####

		# Continue adding more alignments

		# Get the alignment region's kmers
		ref_k_set = getKmers(ref_record.seq[start:end], k_size)

		# score alignments with different weighting schemes
		for sch in schemes:
			x = score(curr_read.k_set, ref_k_set, sch, k_size)
			alignment.scores[sch] = x
		#####
		print("%s\t%s" % (read_name, alignment.toString()))
		# curr_read.alignments.append(alignment)
		# alignments.append(alignment)
	#####


	# For each read, pick the one with the highest score to be primary for each weighting scheme


	# Evaluate true and false rates












if __name__ == "__main__": main()




