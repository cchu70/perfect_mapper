# This script calculates the mapQ based on the same procedure as minimap2 (Heng Li), but uses the number of shared kmers
# mapQsck = 1 - f2/f1
# Where fi is the number of shared unique kemrs between a read and its alignment i. For each read, we only consider the first (f1) and second best (f2)
import sys
import os.path
import math


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
	MQ        = 0.0
	shared_sck_count = 0.0
	order_score = 0.0
	total_shared_sck_count = 0.0
	align_type = None
	data = ""

	# To score based on the unique kmer counts
	score = 0

	def set(self, read_name, start_idx, end_idx, MQ, shared_sck_count, order_score, align_type, data_string, total_shared_sck_count):
		self.read_name = read_name
		self.start_idx = start_idx
		self.end_idx = end_idx
		self.MQ = MQ
		self.shared_sck_count = shared_sck_count
		self.order_score = order_score
		self.align_type = align_type
		self.data = data_string
		self.total_shared_sck_count = total_shared_sck_count
		# print("Made alignment for %s; align_type: %s, total_shared_sck_count: %d, shared_sck_count: %d" % (self.read_name, self.align_type, self.total_shared_sck_count, self.shared_sck_count))
	#####

	def mapQScore(self):
		# print(self.total_shared_sck_count)
		if (self.total_shared_sck_count != 0.0):
			
			return 40 * (self.shared_sck_count / self.total_shared_sck_count) * min(1, self.order_score/500)*self.shared_sck_count
		else:
			return 0.0
		#####
		

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
		return "%s\t%0.5f\t%d\t%d\t%s" % (self.read_name, self.MQ, self.shared_sck_count, self.order_score, self.align_type)

	def print_score(self):
		print("%s\t%0.5f" % (self.data, self.score)) 


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

def which_align_type(string):
	if (string == "tp:A:P"): 
		return "primary"
	elif (string == "tp:A:S"): 
		return "secondary"
	elif (string == "tp:A:I"): 
		return "inversion"
	else:
		print("no alignment type for %s" % string)
		assert False

####

class PAFAlign(AlignData):

	# Somehow allow dynamic selection of these indices
	def __init__(self, paf_string):
		data = paf_string.split()

		read_name = data[0]
		start_idx = int(data[6])
		end_idx = int(data[7])
		MQ = float(data[11])
		shared_sck_count = float(data[-1])
		order_score = float(data[-2])
		total_shared_sck_count = float(data[-3])
		align_type = which_align_type(data[12])

		self.set(read_name, start_idx, end_idx, MQ, shared_sck_count, order_score, align_type, paf_string, total_shared_sck_count)
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
		align_data = PAFAlign(line.strip())

		try:
			alignments[align_data.read_name].append(align_data)
		except KeyError:
			# Assume each line is independent in a PAF file
			# alignments[align_data.read_name] = ReadAlignments(align_data)
			alignments[align_data.read_name] = [align_data]
		#####
	#####

	for read_name in alignments:
		for read_align in alignments[read_name]:
			read_align.score = read_align.mapQScore()
			read_align.print_score()
		#print(read_aligns)
		#print(read_aligns)
		# if(toptwo.second_best != None):
		# 	print(toptwo.second_best)
	#####






if __name__ == "__main__": main()