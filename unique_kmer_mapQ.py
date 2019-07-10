# This script calculates the mapQ based on the same procedure as minimap2 (Heng Li), but uses the number of shared kmers
# mapQsck = 1 - f2/f1
# Where fi is the number of shared unique kemrs between a read and its alignment i. For each read, we only consider the first (f1) and second best (f2)
import sys
import os.path


class TopTwo:

	best = None
	second_best = None

	def add(self, align_data):
		# print(align_data.data)
		if (self.best == None):
			self.best = align_data

		elif (align_data.greaterThan(self.best)):
			# print("Original best %s" % self.best.sck_count)
			self.second_best = self.best
			self.best = align_data
			# print("New best: %d" % align_data.sck_count)
		elif (align_data.greaterThan(self.second_best)):
			self.second_best = align_data
			#print("Second best: %d" % align_data.sck_count)
		#####
	#####

	def __init__(self, align_data_A, align_data_B=None):

		if (align_data_A.greaterThan(align_data_B)):
			self.best = align_data_A
			self.second_best = align_data_B
		else:
			self.best = align_data_B
			self.second_best = align_data_A
		#####
		# print("Best: %s\n sck_count %s" % (self.best.data, self.best.sck_count))\

	def mapQScore(self):
		if (self.second_best != None):
			if (self.best.sck_count > 0):
				self.best.score =  1.0 - (float(self.second_best.sck_count)/float(self.best.sck_count))
			else:
				self.best.score = 0.0
		else:
			self.best.score =  1.0
		#####
#####

class ReadAlignments:
	primary = None # Type alignData
	secondary = None

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
				self.primary.incr_sck_count(align_data.sck_count)
			else:
				self.primary = align_data
			#####
		else:
			if (self.secondary):
				self.secondary.incr_sck_count(align_data.sck_count)
			else:
				self.secondary = align_data
			#####
		#####
	#####

	# def mapQScore(self):
		# 40 * (1 - f2/f1) * min(1, m/200) * log f1



class AlignData:
	read_name = ""
	sck_count = 0
	order_score = 0
	MQ = 0

	data = ""
	score = 0.0

	is_primary = False

	supp = None

	def __init__(self, data_string, sck_count, score=0):
		self.sck_count = sck_count
		self.data = data_string
		self.score = score
		self.primary = self.isPrimary(self.split(data_string)["flag"])
	#####

	def greaterThan(self, align_data):
		if (align_data == None):
			return True
		else:
			return self.sck_count > align_data.sck_count 
		#####
	#####

	def set_sck_count(self,i):
		self.sck_count = i
	#####

	def incr_sck_count(self, to_add):
		self.sck_count = self.sck_count + to_add
	#####

	def __str__(self):
		return "%s\t%0.5f" % (self.data, self.score)
#####


class Minimap2Alignment(AlignData):

	flag = 0
	start_idx = 0
	end_idx = 0

	def split(self, bam_string):
		data = bam_string.split()
		self.read_name = data[0]
		self.flag = int(data[1])
		self.start_idx = int(data[3])
		self.end_idx = int(data[4])
		

	def isPrimary(self,flag):
		if (flag & 256 == 0):
			# primary alignment
			return True
		else:
			# secondary alignment
			return False
		#####
	#####

def main():

	# Get the file with the read, the alignment, original mapping score, and the number of shared unique kmers

	map_sck_counts = sys.argv[1]

	if (not os.path.isfile(map_sck_counts)):
		sys.stderr.write("%s does not exist\n" % map_sck_counts)
	else:
		sys.stderr.write("Loading file %s\n" % map_sck_counts)
	#####
	read_name_idx = int(sys.argv[2])
	sck_count_idx = int(sys.argv[3])

	
	alignments = {}
	# Collect all the shared unique kmers counts

	for line in open(map_sck_counts, "r"):
		# Dictionary of reads maintains queue of the sck alignment scores for each read
		data = line.split()
		read_name = data[read_name_idx]
		sck_count = int(data[sck_count_idx])
		align_data = Minimap2Alignment(line.strip(), float(sck_count))

		try:
			alignments[read_name].add(align_data)
		except KeyError:
			alignments[read_name] = ReadAlignments(align_data)
		#####
	#####

	for read_name in alignments:
		toptwo = alignments[read_name]
		toptwo.mapQScore()
		
		print(toptwo.best)
		# if(toptwo.second_best != None):
		# 	print(toptwo.second_best)
	#####






if __name__ == "__main__": main()