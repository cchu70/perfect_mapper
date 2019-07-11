# Script to calculate the order score based on unique kmers


import sys
import os.path


#=================================================================================
# Functions
#=================================================================================

# Scoring methods
def score_correct_order_only(start, end, read_kmer_idx):
	# Go through each alignment of the previous read
	# sys.stderr.write("In score_correct_order_only\nStart: %d, End: %d\n" %(start, end))
	score = 0
	begin = False
	ref_idx_prev = 0

	for ref_idx in read_kmer_idx:

		if (ref_idx >= start and ref_idx <= end):
			if (not begin):
				# Intitialize
				ref_idx_prev = ref_idx
				begin = True
			else:
				if (ref_idx > ref_idx_prev):
					score = score + 1
				#####

				# Update
				ref_idx_prev = ref_idx
			#####
		#####
	#####
	return score
#####

def score_correct_incorrect_ordering(start, end, read_kmer_idx):
	# Go through each alignment of the previous read
	score = 0
	begin = False
	ref_idx_prev = 0

	for ref_idx in read_kmer_idx:

		if (ref_idx >= start and ref_idx <= end):
			if (not begin):
				# Intitialize
				ref_idx_prev = ref_idx
				begin = True
			else:
				if (ref_idx > ref_idx_prev):
					score = score + 1
				else:
					score = score - penalty
				#####

				# Update
				ref_idx_prev = ref_idx
			#####
		#####
	#####
	return score
#####


# Return the number of shared unique kmers with the aligned region
def count_shared_sck(start, end, read_kmer_idx):
	# sys.stderr.write("In count_shared_sck\n")
	score = 0
	begin = False

	for ref_idx in read_kmer_idx:

		if (ref_idx >= start and ref_idx <= end):
			score = score + 1
			#####
		#####
	#####

	return score


# Other functions

def parseDump(dump_filename):
	kmer_indices = {}	
	with open(dump_filename, "r") as dh:
		curr_read = ""
		ref_indices = []
		for line in dh:
			dump_data   = line.split()
			read_name = dump_data[0]
			ref_idx = dump_data[-1]
			if (not curr_read): 
				# Initialize
				curr_read = read_name
				ref_indices.append(int(ref_idx))
				continue
			elif (curr_read == read_name):
				# Update data
				ref_indices.append(int(ref_idx))
			else:
				# Switching to a new read
				
				kmer_indices[curr_read] = ref_indices
				ref_indices = []
				curr_read = read_name
			#####
		#####
		kmer_indices[curr_read] = ref_indices
	#####
	return kmer_indices

#=================================================================================
# Definitions
#=================================================================================
score_types_func = {'plus1_only': score_correct_order_only,
					'plus_minus': score_correct_incorrect_ordering,
					'mapQ': count_shared_sck}



penalty = 0
reward = 1
args = None

#=================================================================================
# Main
#=================================================================================
def main():
	# Retrieve the dump file (read_name <tab> read_index <tab> ref_index)
	dump = sys.argv[1]

	# Retrieve the alignment file
	# Specify columns later
	maps = sys.argv[2]

	# Sorted and merged bed file of each read's alignments
	align_merge_bed = sys.argv[3]

	if (not os.path.isfile(dump)):
		sys.stderr.write("%s is not a file. Halting execution" % dump)
		assert False
	#####

	if (not os.path.isfile(maps)):
		sys.stderr.write("%s is not a file. Halting execution" % maps)
		assert False
	#####


	idx_read_name = sys.argv[4]
	idx_start = int(sys.argv[5])
	idx_end = int(sys.argv[6])

	score_type = sys.argv[7]

	if (score_type not in score_types_func):
		sys.stderr.write("%s is not a valid scoring option. Select from the following: \n %s\n" % (score_type, '\n\t-'.join(score_types_func)))
		assert False
	else:

		scoreMashMapAlignments = score_types_func[score_type]
	#####

	sys.stderr.write("%s,%s,%s,%s"%(idx_read_name, idx_start, idx_end,score_type))

	# For each read, get the indices of the unique kmers
	kmer_indices = parseDump(dump)

 	# Count the total number of unique kmers on all the alignments for a single read (overlapping alignments)
	total_shared = {}
	with open(align_merge_bed, "r") as bh:
		for line in bh:
			data = line.split()
			read_name = data[0]
			start = float(data[1])
			end = float(data[2])
			try:
				read_kmer_idx = kmer_indices[read_name]
				sck_count = count_shared_sck(start, end, read_kmer_idx)
				total_shared[read_name] = total_shared[read_name] + sck_count # Add all the merged regions shared unique kmers together
			except KeyError:
				read_kmer_idx = kmer_indices[read_name]
				sck_count = count_shared_sck(start, end, read_kmer_idx)
				total_shared[read_name] = sck_count
			#####
		#####
	sys.stderr.write("Done with bed file")

	
	# Get the actual sck counts for each alignment
	with open(maps, "r") as mh:
		#print("# HEADER # Last two columns: %s_score\ttotal_shared_sck" % (score_type))
		for line in mh:
			data = line.split()
			read_name = data[0] 
			start = int(data[idx_start])
			end = int(data[idx_end])
			try:
				read_kmer_idx = kmer_indices[read_name]
				score = scoreMashMapAlignments(start, end, read_kmer_idx)
				print("%s\t%d\t%d" % (line.strip(), score, total_shared[read_name]))
			except KeyError:
				pass
			#####
		#####

		# # Reached end of file
		# try:
		# 	read_kmer_idx = kmer_indices[read_name]
		# 	score = scoreMashMapAlignments(start, end, read_kmer_idx)		
		# 	print("%s\t%d\t%d" % (line.strip(), score, total_shared[read_name]))
		# except KeyError:
		# 	pass
		# #####

	



	

#####




if __name__ == "__main__": main()






	