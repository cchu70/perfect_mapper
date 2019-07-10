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
		if (self.second_best):
			self.best.score =  1.0 - (float(self.second_best.sck_count)/float(self.best.sck_count))
		else:
			self.best.score =  1.0
		#####
#####

class AlignData:
	sck_count = 0
	data = ""
	score = 0.0

	def __init__(self, data_string, sck_count, score=0):
		self.sck_count = sck_count
		self.data = data_string
		self.score = score
	#####

	def greaterThan(self, align_data):
		if (align_data == None):
			return True
		else:
			return self.sck_count > align_data.sck_count 
		#####
	#####

	def __str__(self):
		return "%s\t%0.5f" % (self.data, self.score)
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

	# sys.stderr.write("Read name index: %d\nSCK count index: %d\n" %(read_name_idx, sck_count_idx))

	# For each read, get the top two alignments if there are at least 2, else just take the only one

	alignments = {}

	# Get the top two alignments for each read

	for line in open(map_sck_counts, "r"):
		# Dictionary of reads maintains queue of the sck alignment scores for each read
		data = line.split()
		read_name = data[read_name_idx]
		sck_count = int(data[sck_count_idx])
		align_data = AlignData(line.strip(), float(sck_count))
		#print(align_data)
		# print("Read Name: %s\t sck_count: %d" % (read_name, sck_count))

		try:
			alignments[read_name].add(align_data)
		except KeyError:
			alignments[read_name] = TopTwo(align_data)
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