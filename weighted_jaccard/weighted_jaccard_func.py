# This Script outlines the major functions used in weighted_jaccard_main.py

# Counts the number in each set
def counts(read_k_set, align_k_set, k_size):

	shared_unique_sum = 0
	shared_non_unique_sum = 0

	non_shared_unique_sum = 0
	non_shared_non_unique_sum = 0
	total_sum = 0
	
	# Query the read onto the align set. If match, mark true and increment
	# Sets initialized so the kmers are the keys and all have true as the value

	for k in read_k_set:
		if k in align_k_set:
			if k.isUnique():
				shared__unique_sum += 1
			else:
				shared_non_unique_sum += 1
			#####
			align_k_set[k] = False
			read_k_set[k] = False
		else:
			if k.isUnique():
				non_shared_unique_sum += 1
			else:
				non_shared_non_unique_sum += 1
			#####
		#####
	#####

	for k in align_k_set:
		if align_k_set[k]:
			if k.isUnique():
				non_shared_unique_sum += 1
			else:
				non_shared_non_unique_sum += 1
			#####
		#####
	#####

	return shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum

def weightJaccard(w_unique, w_non_unique, shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum):

	intersection = w_unique * shared_unique_sum + w_non_unique * shared_non_unique_sum
	union = w_unique * non_shared_unique_sum + w_non_unique * non_shared_non_unique_sum

	return float(intersection)/float(union)

def getKmers(seq_str, k_size):
	k_set = {}
	for i in range(len(seq_str) - k_size):
		k = seq_str[i: i + k_size]
		k_set[k] = True
	#####
	return k_set
#####

def score(k_set1, k_set2, sch, k_size):
	shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum = counts(k_set1, k_set2, k_size)
	score = weightJaccard(sch[0], sch[1], shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum)


def parsePaf(paf_string):
	# Based on sam file
	data = paf_string.split()
	read_name = data[0]
	length = int(data[1])
	start = int(data[3])
	cigar = data[-2].split(":")[-1]

	ground_truth = data[-1]

	end = start + parseCigar(cigar)
	

	return read_name, length, start, end, ground_truth

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


###########################################################################
# Classes
###########################################################################

class Read:
	read_name = ""
	length = 0
	k_set = {}
	alignments = []

	def add_alignment(self, alignment):
		self.alignments.append(alignment)
	#####

	def print_alignments(self):
		for alignment in self.alignments:
			print("%s\t%s" % (self.read_name, alignment.toString()))

	def __init__(self, read_name, length, k_set):
		self.read_name = read_name
		self.length = length
		self.k_set = k_set

class Kmer:
	seq = ""
	isUnique = False

class Alignment:
	start_idx = 0
	end_idx = 0
	scores = {}
	ground_truth = False

	def add_score(self, scheme, score):
		self.scores[scheme] = score
	#####

	def toString(self):
		scores_string = "\t".join(["%s:%d" % (sch, score) for sch, score in self.scores])
		return "%d\t%d\t%b\t%s" % (start_idx, end_idx, ground_truth, scores_string)
	#####

	def __init__(self, start, end, ground_truth):
		self.start_idx = start
		self.end_idx = end
		self.ground_truth = ground_truth
#####



