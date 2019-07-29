# Script to procude a true sck counts and ordering score

# Script to merge the supplementary alignments to it's corresponding primary alignment
# Assumes that the sam file has not been altered or sorted. Requires that supplementary alignments immediately follow the main alignment

import sys
import operator


class Read:

	name = ""
	length = 0
	alignments = [] # List of alignments


	total_sck = 0

	def __init__(self, read_name):
		self.read_name = read_name

	def add(self, alignment):
		self.alignments.append(alignment)


	def finish(self):
		for alignment in self.alignments:
			# Calc our mapQ score
			
			score = scoreMapQ(self.total_sck, alignment.shared_sck_count, alignment.order_score)
			toPrint = "%s\t%d\t%d" % (alignment.toString(), self.total_sck, score)
			print(toPrint)
		#####
	#####
			
		
def scoreMapQ(total_sck, shared_sck_count, order_score):
	if (total_sck == 0):
		return 0
	#####
	return 40 * (shared_sck_count / total_sck) * min(1, order_score/500) * shared_sck_count

def sortMerge(alignments):

	regions = [(a.start_idx, a.end_idx) for a in sorted(alignments)]

	tmp_array = []
	x = y = None
	for start, end in regions:
		if not x or not y:
			x = start
			y = end
		else:
			if (start > y):
				tmp_array.append((x,y))
				x = start
				y = end
			else:
				y = end
			#####
		#####
	#####
	tmp_array.append((x,y))
	return tmp_array


class Alignment:

	read_name = ""
	align_type = "S" # I guess I'm making my own alignment file now, either P or S based on minimap
	start_idx = 0
	end_idx = 0
	MQ = 0

	shared_sck_count = 0
	order_score = 0

	def __init__(self, read_name, align_type, start, end, MQ, shared_sck_count, order_score):
		self.read_name = read_name
		self.align_type = align_type
		self.start_idx = start
		self.end_idx = end
		self.MQ = MQ
		self.shared_sck_count = shared_sck_count
		self.order_score = order_score
		

	def __str__(self):
		return "%s\t%s\t%d\t%d\t%d\t%d\t%d" % (self.read_name, self.align_type, self.start_idx, self.end_idx, self.MQ, self.shared_sck_count, self.order_score)

	def toString(self):
		return "%s\t%s\t%d\t%d\t%d\t%d\t%d" % (self.read_name, self.align_type, self.start_idx, self.end_idx, self.MQ, self.shared_sck_count, self.order_score)


	def __cmp__(self, other):
		return self.start_idx - other.start_idx


def parseSam(sam_string):
	# Based on sam file
	read_name = sam_string[0]
	flag =int(sam_string[1])
	start = int(sam_string[3])
	MQ = int(sam_string[4])
	cigar = sam_string[5]

	align_type = alignType(flag)
	align_length = parseCigar(cigar)
	end = start + align_length

	return read_name, align_type, start, end, MQ

def parseCigar(cigar_string):

	length = 0
	num_string = ""

	for c in cigar_string:
		if c.isdigit():
			num_string += c
		else:
			# A letter
			d = int(num_string) # Convert the current num_string into a number
			if (c == "M" or c == "D"):
				# Only add up matches and deletions in the read
				length += d
			#####
			num_string = ""
		#####
	#####
	return length
#####

def alignType(flag):
	if (flag == 0):
		return "P"
	elif (flag & 0x800):
		return "supplementary"
	elif (flag & 256):
		# Secondary alignment
		return "S"
	else:
		# Ignoring reverse alignments
		return "unmapped"
	#####
#####
	
def parseDump(dump_filename):
	kmer_indices = {}	
	with open(dump_filename, "r") as dh:
		curr_read = ""
		visited_idx = {} # track what has been visited
		ref_indices_order = [] # track order
		for line in dh:
			dump_data = line.split()
			read_name = dump_data[0]
			ref_idx = int(dump_data[-1])


			if (not curr_read): 
				# Initialize
				curr_read = read_name
			elif (curr_read != read_name):
				# Switching to a new read
				kmer_indices[curr_read] = (visited_idx, ref_indices_order)
				visited_idx = {}
				ref_indices_order = []
				curr_read = read_name
			#####

			# Update the dictionary
			visited_idx[ref_idx] = 0 # Has not been visited
			ref_indices_order = [ref_idx] # keep the order
		#####
		kmer_indices[curr_read] = visited_idx
	#####
	return kmer_indices


def CountSharedSCKs(alignment_regions, kmer_indices, k_size):
	score = 0
	for start, end in alignment_regions:
		for i in kmer_indices:
			if (i >= start and i < end - k_size):
				score += 1
			#####
		#####
	#####
	return score

def CalcOrderScore(alignment_regions, kmer_indices, k_size):
	score = 0
	prev = sys.maxint # Force skip of the first index
	for start, end in alignment_regions:
		for i in kmer_indices:
			if (i >= start and i < end - k_size and i > prev):
				score += 1
			#####
			prev = i
		#####
	#####
	return score

def score(alignment_regions, kmer_indices, k_size):
	sck_count = 0
	order_score = 0
	prev = sys.maxint # Force skip of the first index
	for start, end in alignment_regions:
		for i in kmer_indices:
			if (i >= start and i < end - k_size):
				sck_count += 1
				if (i > prev):
					order_score += 1
				#####
			#####
			prev = i
		#####
	#####
	return sck_count, order_score

def main():

	# Get sam file
	sam_fh = sys.stdin

	k_size = 21


	# Dump File
	dump = sys.argv[1]

	# For each read, get the indices of the unique kmers
	sys.stderr.write("Parsing dump file %s ...\n" % dump)
	kmer_indices = parseDump(dump)
	sys.stderr.write("Done parsing dump file %s ...\n" % dump)

	# maps = sys.argv[2]

	# initialize starting variables
	curr_align = None
	curr_read = None

	curr_read_alignments = []

	# sys.stderr.write("Reading map file %s\n" % maps)
	# Read through the sam file
	for line in sam_fh:
		read_name, align_type, start, end, MQ = parseSam(line.split())

		# Initialize
		if not curr_read:
			curr_read = Read(read_name)
			sys.stderr.write("Parsing alignments for read %s\n" % read_name)
		#####

		if (align_type == "P" or align_type == "S"):
			# Add the representative alignment

			try:
				shared, order = score([(start, end)], kmer_indices[read_name], 21)
				curr_align = Alignment(read_name, align_type, start, end, MQ, shared, order)

				if (read_name == curr_read.read_name):
					curr_read.add(curr_align)
				else:
					# Calc and print data
					regions = sortMerge(curr_read.alignments)

					curr_read.total_sck = CountSharedSCKs(regions, kmer_indices[curr_read.read_name], k_size)
					curr_read.finish()
					sys.stderr.write("Finished parsing alignments for read %s\n" % read_name)

					# start new read
					sys.stderr.write("Parsing alignments for read %s\n" % read_name)
					curr_read = Read(read_name)
					curr_read.add(curr_align)
			except KeyError:
				pass

		#####

			
	######

	# Finish last one
	# Calc and print data


	regions = sortMerge(curr_read.alignments)

	#print(curr_read.alignments)
	curr_read.total_sck = CountSharedSCKs(regions, kmer_indices[curr_read.read_name], k_size)
	print(curr_read.total_sck)
	curr_read.finish()
	sys.stderr.write("Finished parsing alignments for read %s\n" % read_name)
	

#####


if __name__ == "__main__": main()

