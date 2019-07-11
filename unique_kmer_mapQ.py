# This script calculates the mapQ based on the same procedure as minimap2 (Heng Li), but uses the number of shared kmers
# mapQsck = 1 - f2/f1
# Where fi is the number of shared unique kemrs between a read and its alignment i. For each read, we only consider the first (f1) and second best (f2)
import sys
import os.path
import math
from enum import Enum


#class TopTwo:

	# best = None
	# second_best = None

	# def add(self, align_data):
	# 	# print(align_data.data)
	# 	if (self.best == None):
	# 		self.best = align_data

	# 	elif (align_data.greaterThan(self.best)):
	# 		# print("Original best %s" % self.best.sck_count)
	# 		self.second_best = self.best
	# 		self.best = align_data
	# 		# print("New best: %d" % align_data.sck_count)
	# 	elif (align_data.greaterThan(self.second_best)):
	# 		self.second_best = align_data
	# 		#print("Second best: %d" % align_data.sck_count)
	# 	#####
	# #####

	# def __init__(self, align_data_A, align_data_B=None):

	# 	if (align_data_A.greaterThan(align_data_B)):
	# 		self.best = align_data_A
	# 		self.second_best = align_data_B
	# 	else:
	# 		self.best = align_data_B
	# 		self.second_best = align_data_A
	# 	#####
	# 	# print("Best: %s\n sck_count %s" % (self.best.data, self.best.sck_count))\

	# def mapQScore(self):
	# 	if (self.second_best != None):
	# 		if (self.best.sck_count > 0):
	# 			self.best.score =  1.0 - (float(self.second_best.sck_count)/float(self.best.sck_count))
	# 		else:
	# 			self.best.score = 0.0
	# 	else:
	# 		self.best.score =  1.0
	# 	#####
#####

class ReadAlignments:
	primary = None # Type alignData
	secondary = None
	score = 0

	def __init__(self, prim, sec):
		self.primary = primary
		self.secondary = sec
	#####

	def __init__(self, align_data):
		if (align_data.is_primary):
			self.primary = align_data
		else:
			self.secondary = align_data
		#####
	#####

	def add_align(self, align_data):
		if (align_data.is_primary):
			if (self.primary):
				self.primary.incr_shared_sck_count(align_data.shared_sck_count)
			else:
				self.primary = align_data
			#####
		else:
			if (self.secondary):
				self.secondary.incr_shared_sck_count(align_data.shared_sck_count)
			else:
				self.secondary = align_data
			#####
		#####
	#####

	def mapQScore(self):
		# 40 * (1 - f2/f1) * min(1, m/200) * log f1

		if (self.primary.shared_sck_count == 0):
			self.score =  0
		elif (self.secondary):
			self.score =  40 * (1 - self.secondary.shared_sck_count / self.primary.shared_sck_count) * min(1.0, self.primary.order_score / 200)
			
			print("Secondary exists")
			print(self)
			print(self.score)
		else:
			self.score =  40 *  min (1.0, self.primary.order_score / 200)
			#print(self.primary)
		#####
	#####

	def __str__(self):
		if (self.secondary):
			return "Score: %0.5f\nPrimary: %s\n Secondary: %s" % (self.score, self.primary, self.secondary)
		else:
			return "Score: %0.5f\nPrimary: %s" % (self.score, self.primary)



class AlignData:
	read_name = ""
	start_idx = 0
	end_idx   = 0 
	MQ        = 0
	shared_sck_count = 0
	order_score = 0
	align_type = None
	data = ""

	# To score based on the unique kmer counts
	score = 0

	def greaterThan(self, align_data):
		if (align_data == None):
			return True
		else:
			return self.shared_sck_count > align_data.shared_sck_count 
		#####
	#####

	def set_shared_sck_count(self,i):
		self.shared_sck_count = i
	#####

	def incr_shared_sck_count(self, to_add):
		self.shared_sck_count = self.shared_sck_count + to_add
	#####

	def __str__(self):
		return "%s\t%0.5f" % (self.data, self.shared_sck_count)
#####


class Minimap2Alignment(AlignData):

	flag = 0
	start_idx = 0
	end_idx = 0

	def split(self, bam_string):
		data = bam_string.split()
		self.read_name = data[0]
		self.flag = int(data[1])
		self.is_primary = self.isPrimary(self.flag)

		self.start_idx = float(data[3])
		self.end_idx = float(data[4])
		self.MQ = float(data[5])

		self.shared_sck_count = float(data[6])
		# self.order_score = float(data[3])
		self.order_score = float(data[8])
		self.data = bam_string

		

	def isPrimary(self,flag):
		if (flag & 256 == 0):
			# primary alignment
			return True
		else:
			# secondary alignment
			return False
		#####
	#####

class AlignType(Enum):
	PRIMARY, SECONDARY

class PAFAlign(AlignData):
	align_type = None

	def __init__(self, paf_string):
		data = paf_string.split()

		read_name = data[0]
		start_idx = data[0]
		end_idx = data[0]
		MQ = data[0]
		shared_sck_count = data[0]
		order_score = data[0]
		align_type = data[0]

		self(read_name, start_idx, end_idx, MQ, shared_sck_count, order_score, align_type, paf_string)
	#####

	def __init__(self, read_name, start_idx, end_idx, MQ, shared_sck_count, order_score, align_type, data_string):
		self.read_name = read_name
		self.start_idx = start_idx
		self.end_idx = end_idx
		self.MQ = MQ
		self.shared_sck_count = shared_sck_count
		self.order_score = order_score
		self.align_type = align_type
		self.data = data_string
	#####


def main():

	# Get the file with the read, the alignment, original mapping score, and the number of shared unique kmers

	map_shared_sck_counts = sys.argv[1]

	if (not os.path.isfile(map_shared_sck_counts)):
		sys.stderr.write("%s does not exist\n" % map_shared_sck_counts)
	else:
		sys.stderr.write("Loading file %s\n" % map_shared_sck_counts)
	#####
	# read_name_idx = int(sys.argv[2])
	# shared_sck_count_idx = int(sys.argv[3])
	# order_score_idx = int(sys.argv[4])

	
	alignments = {}
	# Collect all the shared unique kmers counts

	for line in open(map_shared_sck_counts, "r"):
		# Dictionary of reads maintains queue of the sck alignment scores for each read
		# data = line.split()
		# read_name = data[read_name_idx]
		# shared_sck_count = int(data[shared_sck_count_idx])
		# order_score = 
		align_data = Minimap2Alignment(line.strip())

		try:
			alignments[align_data.read_name].add_align(align_data)
		except KeyError:
			alignments[align_data.read_name] = ReadAlignments(align_data)
		#####
	#####

	for read_name in alignments:
		read_aligns = alignments[read_name]
		read_aligns.mapQScore()
		#print(read_aligns)
		#print(read_aligns)
		# if(toptwo.second_best != None):
		# 	print(toptwo.second_best)
	#####






if __name__ == "__main__": main()