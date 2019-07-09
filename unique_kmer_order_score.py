# Script to calculate the order score based on unique kmers


import sys
import os.path


def main():
	# Retrieve the dump file (read_name <tab> read_index <tab> ref_index)
	dump = sys.argv[1]

	# Retrieve the alignment file
	# Specify columns later
	maps = sys.argv[2]

	if (not os.path.isfile(dump)):
		sys.stderr.write("%s is not a file. Halting execution" % dump)
	#####

	if (not os.path.isfile(maps)):
		sys.stderr.write("%s is not a file. Halting execution" % maps)
	#####


	idx_read_name = sys.argv[3]
	idx_start = int(sys.argv[4])
	idx_end = int(sys.argv[5])

	# For each read, get the indices of the unique kmers
	kmer_indices = parseDump(dump)

	with open(maps, "r") as mh:

		for line in mh:
			data = line.split()
			read_name = data[0] 
			try:
				read_kmer_idx = kmer_indices[read_name]
				start = int(data[idx_start])
				end = int(data[idx_end])
				score = scoreMashMapAlignments(read_name, start, end, read_kmer_idx)
				print("%s\t%d\t%d" % (line.strip(), score, len(read_kmer_idx)))
			except:
				continue
		#####

		# Reached end of file
		try:
			read_kmer_idx = kmer_indices[read_name]
			score = scoreMashMapAlignments(read_name, start, end, read_kmer_idx)		
			print("%s\t%d\t%d" % (line.strip(), score, len(read_kmer_idx)))
		except:
			pass
#####

def scoreMashMapAlignments(read_name, start, end, read_kmer_idx):
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
				#####

				# Update
				ref_idx_prev = ref_idx
			#####
		#####
	#####
	return score
	
#####

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



if __name__ == "__main__": main()






	