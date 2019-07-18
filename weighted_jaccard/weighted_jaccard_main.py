# This is the main script to test different weighting schemes to compute a weighted jaccard score on shared kmers between a read and it's alignment region

from  weighted_jaccard_func import getKmers, score, parsePaf, Read, Alignment, parseUniqueFile, parseSam, parseFasta
#from Bio import SeqIO
import sys
import os
import subprocess

def main():

	sch_start = 1
	sch_end = 16
	k_size = 21
	read_fasta = "/data/Phillippy/projects/perfect-polish/chr22_info/chr22.sim_reads.fasta"
	ref_fasta = "/data/Phillippy/projects/perfect-polish/chr22_info/chr22.fasta"
	align_file = "/data/Phillippy/projects/perfect-polish/chr22_info/representative_only.multiple_aligns_only.rev_false.minimap2_N50_30kb.real.sam"
	unique_k_file = "/data/Phillippy/projects/perfect-polish/chr22_info/chr22.asm.sck_list.txt"


	try:
		read_fasta = sys.argv[1]
		align_file = sys.argv[2]
	except:
		pass
	#####

	# Set up schemes
	w = 2
	schemes = [(1, w ** i) for i in range(sch_start, sch_end)]

	# Get read sequences
	sys.stderr.write("Parsing Read fasta: %s\n" % read_fasta)
	read_records = parseFasta(read_fasta) # Dictionary of read names and it's corresponding sequence

	sys.stderr.write("Parsing Ref fasta: %s\n" % ref_fasta)
	ref_record = list(parseFasta(ref_fasta).values())[0] # Should only be one reference

	sys.stderr.write("Parsing unique file: %s\n" % unique_k_file)
	unique_table = parseUniqueFile(unique_k_file)
	sys.stderr.write("Number of unique kmers: %d\n" % len(unique_table))

	curr_read = None

	for line in open(align_file, "r"):

		read_name, length, ref_start, ref_end, ground_truth, read_start, read_end = parseSam(line.strip())

		sys.stderr.write("Length: %d, ref_start: %d, ref_end: %d, truth: %s, read_start: %d, read_end: %d\n" % (length, ref_start, ref_end, ground_truth, read_start, read_end))

		alignment = Alignment(ref_start, ref_end, ground_truth)


		# Check which read (current or next) this alignment corresponds to 
		if (curr_read):
			if (read_name != curr_read.read_name):
				# evaluate the curr read performance
				read_records.pop(curr_read, None) # remove from the sequence dictionary because I am done with it? will this run faster?
				curr_read = Read(read_name, length, read_records[read_name])
			#####
		else:
			# initialize
			curr_read = Read(read_name, length, read_records[read_name])
		#####

		# Continue adding more alignments

		# Get the alignment region's kmers
		ref_k_set = getKmers(ref_record[ref_start:ref_end], k_size)

		# score alignments with different weighting schemes
		for sch in schemes:
			x = score(curr_read.seq_str[read_start:read_end], ref_k_set, sch, k_size, unique_table)
			alignment.scores[sch] = x
		#####
		print("%s\t%s" % (read_name, alignment.toString()))
	#####



if __name__ == "__main__": main()




