# This Script outlines the major functions used in weighted_jaccard_main.py
import time
import sys


def parseUniqueFile(unique_k_file):
	table = {}
	for k in open(unique_k_file, "r"):
		table[k.strip()] = True
	#####
	return table
#####



# Counts the number in each set
def counts(read_seq, align_k_set, k_size, unique_table):

	shared_unique_sum = 0
	shared_non_unique_sum = 0

	non_shared_unique_sum = 0
	non_shared_non_unique_sum = 0
	total_sum = 0

	# Debugging
	k_count = 0
	k_uniq_count = 0
	
	# Query the read onto the align set. If match, mark true and increment
	# Sets initialized so the kmers are the keys and all have true as the value
	for i in range(len(read_seq) - k_size):
		# get read's kmers
		k = str(read_seq[i: i + k_size])

		# Debugging
		k_count += 1

		# sys.stderr.write("%s\n" % isUnique(k))
		if (k not in unique_table):
			print(k)

		if k in align_k_set:
			if k in unique_table:
				# Debugging
				k_uniq_count += 1

				shared_unique_sum += 1
			else:
				shared_non_unique_sum += 1
			#####
			align_k_set[k] = False
		else:
			if k in unique_table:

				# Debugging
				k_uniq_count += 1

				non_shared_unique_sum += 1
			else:
				non_shared_non_unique_sum += 1
			#####
		#####
	#####

	sys.stderr.write("read k_count: %d\n" % k_count)

	for k in align_k_set:

		# Debugging
		k_count += 1
		if align_k_set[k]:
			if k in unique_table:
				non_shared_unique_sum += 1

				# Debugging
				k_uniq_count += 1

			else:
				non_shared_non_unique_sum += 1
			#####
		#####
	#####

	sys.stderr.write("unique k_count union found: %d\n" % k_uniq_count)

	sys.stderr.write("union k_count: %d\n" % k_count)


	return shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum

def weightJaccard(w_non_unique, w_unique, shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum):

	intersection = w_unique * shared_unique_sum + w_non_unique * shared_non_unique_sum
	union = w_unique * non_shared_unique_sum + w_non_unique * non_shared_non_unique_sum + intersection

	return float(intersection)/float(union)

def getKmers(seq_str, k_size):

	k_set = {}
	for i in range(len(seq_str) - k_size):
		k = str(seq_str[i: i + k_size])
		k_set[k] = True
	#####
	return k_set
#####

# def isUnique(k_str):
# 	try:
# 		return unique_table[k_str] 
# 	except KeyError:
# 		return False
# 	#####

def score(read_seq, align_k_set, sch, k_size, unique_table):
	# Debugging
	sys.stderr.write("Non unique weight: %s, unique weight: %s\n" % (sch[0], sch[1]))

	shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum = counts(read_seq, align_k_set, k_size, unique_table)

	# Debugging
	sys.stderr.write("Shared unique: %d, shared non-unique: %d, non-shared unique: %d, non shared non unique: %d\n" % (shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum))
	
	non_unique_weight = sch[0]
	unique_weight = sch[1]


	similarity_score = weightJaccard(non_unique_weight, unique_weight, shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum)
	
	sys.stderr.write("score: %0.5f\n" % (similarity_score))

	assert False
	return score

def parsePaf(paf_string):
	# Based on paf file
	data = paf_string.split()
	read_name = data[0]
	length = int(data[1])
	start = int(data[7])
	end = int(data[8])

	ground_truth = data[-1]

	return read_name, length, start, end, ground_truth

def parseSam(sam_str):
	data = sam_str.split()

	read_name = data[0]

	start = int(data[3])
	cigar = data[5]



	length = parseCigar(cigar)

	end = start + length

	ground_truth = data[-1]

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

			if (c == "M" or c == "D" or c == "N"):
				# Only add up matches and deletions in the read
				length += d
			#####
			num_string = ""
		#####
	#####
	return length
#####

def parseFasta(fasta_file):
	read_record = {}
	header = None
	seq_str = ""
	for line in open(fasta_file, "r"):
		if ">" in line:

			if header:
				read_record[header] = seq_str
				seq_str = ""
			#####

			# Start next header
			header = line.strip().replace(">", "")
		else:
			seq_str += line.strip()
		#####
	#####
	read_record[header] = seq_str
	return read_record



###########################################################################
# Classes
###########################################################################

class Read:
	read_name = ""
	length = 0
	seq_str = ""
	alignments = []

	def add_alignment(self, alignment):
		self.alignments.append(alignment)
	#####

	def print_alignments(self):
		for alignment in self.alignments:
			print("%s\t%s" % (self.read_name, alignment.toString()))

	def __init__(self, read_name, length, seq_str):
		self.read_name = read_name
		self.length = length
		self.seq_str = seq_str


class Alignment:
	start_idx = 0
	end_idx = 0
	scores = {}
	ground_truth = False

	def add_score(self, scheme, score):
		self.scores[scheme] = score
	#####

	def toString(self):
		scores_string = "\t".join(["%s=%0.5f" % (sch, self.scores[sch]) for sch in self.scores])
		return "%d\t%d\t%s\t%s" % (self.start_idx, self.end_idx, self.ground_truth, scores_string)
	#####

	def __init__(self, start, end, ground_truth):
		self.start_idx = start
		self.end_idx = end
		self.ground_truth = ground_truth
#####

class Scheme:
	unique_weight = 1
	non_unique_weight = 1

	def __init__(self, w1, w2):
		self.unique_weight = w1
		self.non_unique_weight = w2
	#####
#####




