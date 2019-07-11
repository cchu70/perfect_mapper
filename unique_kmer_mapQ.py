# This script calculates the mapQ based on the same procedure as minimap2 (Heng Li), but uses the number of shared kmers
# mapQsck = 1 - f2/f1
# Where fi is the number of shared unique kemrs between a read and its alignment i. For each read, we only consider the first (f1) and second best (f2)
import sys
import os.path
import math

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
	length = 0
	start_idx = 0
	end_idx   = 0 
	MQ        = 0.0
	shared_sck_count = 0.0
	order_score = 0.0
	total_shared_sck_count = 0.0
	align_type = None
	data = ""
	isTrue = False

	# To score based on the unique kmer counts
	score = 0

	def set(self, read_name, length, start_idx, end_idx, MQ, shared_sck_count, order_score, align_type, data_string, total_shared_sck_count):
		self.read_name = read_name
		self.length = length
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
			return self.score > align_data.score 
		#####
	#####

	def __str__(self):
		return "%s\t%0.5f\t%d\t%d\t%s" % (self.read_name, self.MQ, self.shared_sck_count, self.order_score, self.align_type)

	def print_score(self):
		return "%s\t%0.5f" % (self.data, self.score)


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
	isPrimary = False

	# Somehow allow dynamic selection of these indices
	def __init__(self, paf_string):
		data = paf_string.split()

		read_name = data[0]
		length = int(data[1])
		start_idx = int(data[7])
		end_idx = int(data[8])
		MQ = float(data[11])
		shared_sck_count = float(data[-1])
		order_score = float(data[-2])
		total_shared_sck_count = float(data[-3])
		align_type = which_align_type(data[12])

		if (align_type == "primary"):
			self.isPrimary = True
		#####

		self.set(read_name, length, start_idx, end_idx, MQ, shared_sck_count, order_score, align_type, paf_string, total_shared_sck_count)
	#####

def parse_bed(bed_file):
	true_origins = {}
	for line in open(bed_file, "r"):
		read_name, start, end = line.strip().split()
		true_origins[read_name] = int(start)
	#####
	return true_origins


def main():

	# Get the file with the read, the alignment, original mapping score, and the number of shared unique kmers

	map_shared_sck_counts = sys.argv[1]

	if (not os.path.isfile(map_shared_sck_counts)):
		sys.stderr.write("%s does not exist\n" % map_shared_sck_counts)
	else:
		sys.stderr.write("Loading file %s\n" % map_shared_sck_counts)
	#####
	
	mapQ_output = sys.argv[2]
	fh = open(mapQ_output, "w")

	
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


	# Collected the primary read according to our calculated mapQ score
	mapQBest = {}

	for read_name in alignments:
		best_align = None
		for read_align in alignments[read_name]:
			read_align.score = read_align.mapQScore()
			fh.write(read_align.print_score())

			if (read_align.greaterThan(best_align)):
				best_align = read_align
			#####
		#####

		# Record the best alignment for each read based on our score
		mapQBest[read_name] = best_align
	#####

	fh.close()

##############################################################################
	# Verify correctness
	# Pass in the read's true region

	read_org_bed = sys.argv[3]

	true_origins_start = parse_bed(read_org_bed)

	# For each alignment, check if the alignment is true, and if it is primary or not
	minimap_prim_true_count = 0
	sck_mapQ_prim_true_count = 0
	minimap_prim_false_count = 0
	sck_mapQ_prim_false_count = 0

	minimap_sec_true_count = 0
	sck_mapQ_sec_true_count = 0

	minimap_sec_false_count = 0
	sck_mapQ_sec_false_count = 0

	minimap_sck_mapQ_agree_true = 0
	minimap_sck_mapQ_agree_false = 0

	for read_name in alignments:
		true_start = true_origins_start[read_name]
		best_align = mapQBest[read_name]
		# print(read_name)
		for read_align in alignments[read_name]:
			# 50% covers
			l_bound = true_start - read_align.length / 2.0
			u_bound = true_start + read_align.length / 2.0

			# print(read_align.length)
			# print(true_start)
			# print(l_bound)
			# print(u_bound)
			# print(read_align.start_idx)
			# print("============")
			
			if l_bound <= read_align.start_idx <= u_bound:
				read_align.isTrue = True

				if (read_align.isPrimary and (read_align == best_align)):
					minimap_prim_true_count += 1
					sck_mapQ_prim_true_count += 1
					minimap_sck_mapQ_agree_true += 1
				elif(read_align.isPrimary):
					minimap_prim_true_count += 1
					sck_mapQ_sec_true_count += 1
				elif(read_align == best_align):
					sck_mapQ_prim_true_count += 1
					minimap_sec_true_count += 1
				else:
					# Secondary read for both minimap and our scoring system
					sck_mapQ_sec_true_count += 1
					minimap_sec_true_count += 1
				#####

			else:
				if (read_align.isPrimary and (read_align == best_align)):
					minimap_prim_false_count += 1
					sck_mapQ_prim_false_count += 1
					minimap_sck_mapQ_agree_false += 1
				elif(read_align.isPrimary):
					minimap_prim_false_count += 1
					sck_mapQ_sec_false_count += 1
				elif(read_align == best_align):
					sck_mapQ_prim_false_count += 1
					minimap_sec_false_count += 1
				else:
					# Both methods indicate false secondary alignment
					minimap_sec_false_count += 1
					sck_mapQ_sec_false_count += 1
				#####
			#####
		#####
	###

	print("=======================SUMMARY=======================")
	print("Minimap Primary True: %d" % minimap_prim_true_count)
	print("Minimap Primary False: %d" % minimap_prim_false_count)
	print("Minimap Secondary True: %d" % minimap_sec_true_count)
	print("Minimap Secondary False: %d\n" % minimap_sec_false_count)

	print("MapQ Primary True: %d" % sck_mapQ_prim_true_count)
	print("MapQ Primary False: %d" % sck_mapQ_prim_false_count)
	print("MapQ Secondary True: %d" % sck_mapQ_sec_true_count)
	print("MapQ Secondary False: %d\n" % sck_mapQ_sec_false_count)
	
	print("Minimap and MapQ Primary True: %d" % minimap_sck_mapQ_agree_true)
	print("Minimap and MapQ Primary False: %d" % minimap_sck_mapQ_agree_false)

	print("=====================================================")








if __name__ == "__main__": main()