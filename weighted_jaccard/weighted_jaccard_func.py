# This Script outlines the major functions used in weighted_jaccard_main.py
import time
import sys
from optparse import OptionParser

stages_list = ['start', 'end']


def init_parser():
	parser = OptionParser(usage="usage: %prog [options] filename",
                          version="%prog 1.0")

	parser.add_option("-s", "--stage", default="start", help="Indicate which stage you want to start at in the pipeline", choices=stages_list, dest="start_stage")
	parser.add_option("-c", "--config", help="Input path to configuration file", dest="config_file")
	return parser


# def parseUniqueFile(unique_k_file):
# 	table = {}
# 	for k in open(unique_k_file, "r"):
# 		table[k.strip()] = True
# 	#####
# 	return table
# #####

def parseKmerFile(unique_k_file):
	table = {}
	for line in open(unique_k_file, "r"):
		k, value = line.strip().split()
		table[k] = int(value)
	#####
	return table
#####


# Counts the number in each set
def counts(read_k_set, align_k_set, kmer_table):

	shared_unique_sum = 0
	shared_non_unique_sum = 0
	shared_error_sum = 0

	non_shared_unique_sum = 0
	non_shared_non_unique_sum = 0
	non_shared_error_sum = 0

	total_sum = 0
	
	# Query the read onto the align set. If match, mark true and increment
	# Sets initialized so the kmers are the keys and all have true as the value
	for k in read_k_set:

		# # sys.stderr.write("%s\n" % isUnique(k))
		# if (k not in unique_table):
		# 	print(k)

		if k in align_k_set:
			if k in kmer_table:
				if (kmer_table[k] == 1):
					shared_unique_sum += 1
				else:
					shared_non_unique_sum += 1
			else:
				shared_error_sum += 1
			#####
			# Don't consider this kmer again
			align_k_set[k] = False
		else:
			if k in kmer_table:

				if (kmer_table[k] == 1):
					non_shared_unique_sum += 1
				else:
					non_shared_non_unique_sum += 1
				#####

			else:
				non_shared_error_sum += 1
			#####
		#####
	#####


	# Remaining kmers in the alignment
	for k in align_k_set:
		if align_k_set[k]:
			if k in kmer_table:

				if (kmer_table[k] == 1):
					non_shared_unique_sum += 1
				else:
					non_shared_non_unique_sum += 1
				#####

			else:
				non_shared_error_sum += 1
			#####
		#####
	#####

	return shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum, shared_error_sum, non_shared_error_sum

def weightJaccard(w_non_unique, w_unique, shared_unique_sum, shared_non_unique_sum, non_shared_unique_sum, non_shared_non_unique_sum, shared_error_sum, non_shared_error_sum):

	
	intersection = w_unique * shared_unique_sum + w_non_unique * (shared_non_unique_sum + shared_error_sum)

	
	union = w_unique * non_shared_unique_sum + w_non_unique * (non_shared_non_unique_sum + non_shared_error_sum) + intersection

	# equation = "(%0.3f*%d+%0.3f*%d)/(%0.3f*%d+%0.3f*%d+%d)" % (w_unique, shared_unique_sum, w_non_unique, shared_non_unique_sum, w_unique, non_shared_unique_sum, w_non_unique, non_shared_non_unique_sum, intersection)

	# return (float(intersection)/float(union), equation)
	try:
		x = float(intersection)/float(union)
	except ZeroDivisionError:
		sys.stderr.write("Itersection = %d * %d + %d * %d\n" % (w_unique, shared_unique_sum, w_non_unique, shared_non_unique_sum))
		sys.stderr.write("Union = %d * %d + %d * %d + %d" % (w_unique, non_shared_unique_sum, w_non_unique, non_shared_non_unique_sum, intersection))
		sys.stderr.write("%d\t%d\n" % (intersection, union))
		return -1.0
	#####

	return float(intersection)/float(union)

def getKmers(seq_str, k_size):

	k_set = {}
	for i in range(len(seq_str) - k_size):
		k = str(seq_str[i: i + k_size])
		k_set[k] = True
	#####
	return k_set
#####

def parsePaf(paf_string):
	# Based on paf file
	data = paf_string.split()
	read_name = data[0]
	length = int(data[1])
	start = int(data[7])
	end = int(data[8])

	align_length = int(data[10])

	ground_truth = data[-1]

	return read_name, length, start, end, ground_truth

def parseSam(sam_str):
	data = sam_str.split()

	read_name = data[0]

	ref_start = int(data[3])
	cigar = data[5]

	length, read_start, read_end = parseCigar(cigar)

	if length < 0:
		sys.stderr.write("Unable to parse cigar string %s. Original string was: \n%s\n" % (cigar, sam_str))	
		assert False

	ref_end = ref_start + length

	ground_truth = data[-1]

	map_truth = "P" if int(data[1]) & 256 == 0 else "S"

	return read_name, length, ref_start, ref_end, ground_truth, read_start, read_end, map_truth


def parseMashMap(mashmap_str):
	data = mashmap_str.split()
	read_name = data[0]
	length = int(data[1])
	read_start = int(data[2])
	read_end = int(data[3])
	ref_start = int(data[7])
	ref_end = int(data[8])
	ground_truth = data[-1]
	map_truth = data[-2]
	
	return read_name, length, ref_start, ref_end, ground_truth, read_start, read_end, map_truth

def parseCigar(cigar_string):

	length = 0
	num_string = ""

	read_start = 0
	read_end = None

	for c in cigar_string:
		if c.isdigit():
			num_string += c
		else:
			# A letter
			try:
				d = int(num_string) # Convert the current num_string into a number

				if (c == "M" or c == "D" or c == "N"):
					# Only add up matches and deletions in the read
					length += d

				elif (c == "S" or c == "H"):
					# Get the start and end of the read
					if not read_start:
						read_start = d
					else:
						read_end = -d # Get the index from the back of the read
					#####
				#####
			except:
				sys.stderr.write("Error parsing cigar string : %s\n" % cigar_string)
				return -1, 0, 0
			#####


			num_string = ""
		#####
	#####
	return length, read_start, read_end
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

def parseWJ(wj_str):
	data = wj_str.split("\t")
	read_name = data[0]
	map_truth = data[1]
	start = int(data[2])
	end = int(data[3])
	ground_truth = data[4]
	pid = float(data[5])

	scores = data[6:]
	scores_table = {}
	for s in scores:
		# scheme, score, equation = s.split("=")
		# scores_table[scheme] = (float(score), equation)
		scores_table[s.split("=")[0]] = float(s.split("=")[1])

	return read_name, map_truth, start, end, ground_truth, pid, scores_table

# Constants

align_file_parser = {"sam": parseSam, 
					"mashmap": parseMashMap,
					"paf": parsePaf}

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
	read_name = ""
	map_truth = "S"
	start_idx = 0
	end_idx = 0
	scores = {}
	ground_truth = False
	pid = 0
	data = "" # original string

	def add_score(self, scheme, score):
		self.scores[scheme] = score
	#####

	def scoreString(self):
		scores_string = "\t".join(["(%0.3f, %0.3f)=%s" % (sch.non_unique_weight, sch.unique_weight, self.scores[sch]) for sch in self.scores])
		return "%s\t%s\t%d\t%d\t%s\t%0.8f\t%s" % (self.read_name, self.map_truth, self.start_idx, self.end_idx, self.ground_truth,  self.pid, scores_string)
	#####

	def __init__(self, read_name, map_truth, start, end, ground_truth, pid):
		self.read_name = read_name
		self.map_truth = map_truth
		self.start_idx = start
		self.end_idx = end
		self.ground_truth = ground_truth
		self.pid = pid

	def __str__(self):
		return "%s\t%s\t%d\t%d\t%s\t%0.8f" % (self.read_name, self.map_truth, self.start_idx, self.end_idx, self.ground_truth, self.pid)
#####

class Scheme:
	unique_weight = 1
	non_unique_weight = 1

	def __init__(self, w1, w2):
		self.non_unique_weight = w1
		self.unique_weight = w2
	#####

	def __str__(self):
		return "(%0.5f,%0.5f)" % (self.non_unique_weight, self.unique_weight)
#####




