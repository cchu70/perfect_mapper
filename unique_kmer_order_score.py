# Script to calculate the order score based on unique kmers


import sys


def main():
	# Retrieve the dump file (read_name <tab> read_index <tab> ref_index)
	dump = sys.argv[1]

	# Retrieve the alignment file
	# cat chr22.k21_s500_none.true.out | awk '{print $1"\t"$8"\t"$9}'
	mashmap = sys.argv[2]


	# For each read, find the corresponding alignments *(array of tuples)
	alignments = parseMashMap(mashmap)

	with open(dump, "r") as dh:

		curr_read = ""
		ref_indices = []

		for line in dh:
			dump_data   = line.split()
			read_name = dump_data[0]
			ref_idx = dump_data[-1]
			if (not curr_read): 
				# Initialize
				curr_read = read_name
				continue
			elif (curr_read == read_name):
				# Update data
				ref_indices.append(int(ref_idx))
			else:
				# Switching to a new read

				# Get the scores of the current read
				try:
					read_alignments = alignments[curr_read]
					scoreMashMapAlignments(ref_indices, read_name, read_alignments)
				except:
					sys.stderr.write("%s has not mashmap mapping\n" % curr_read)
				#####

				# Begin the new read
				ref_indices = []
				curr_read = read_name
			#####
		#####

		# Reached end of file
		scoreMashMapAlignments(ref_indices, read_name, alignments[read_name])		
#####

def parseMashMap(mashmap_filename):
	alignments = {}
	with open(mashmap_filename, "r") as mh:
		for line in mh:
			data = line.split()
			read_name = data[0] 
			start = int(data[7])
			end = int(data[8])
			aRange = (start, end, line.strip())
			try: 
				alignments[read_name].append(aRange)
			except:
				alignments[read_name] = [aRange]
			#####
		#####
	#####
	return alignments
#####

def scoreMashMapAlignments(ref_indices, read_name, read_alignments):
	# Go through each alignment of the previous read
	for start, end, mashmap_info in read_alignments:
		score = 0
		ref_idx_prev= 0
		begin = False

		for ref_idx in ref_indices:

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
		print("%s\t%d" % (mashmap_info, score))
#####


if __name__ == "__main__": main()






	