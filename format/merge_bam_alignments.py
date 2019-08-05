# Script to merge the supplementary alignments to it's corresponding primary alignment
# Assumes that the sam file has not been altered or sorted. Requires that supplementary alignments immediately follow the main alignment

import sys



class Alignment:

	read_name = ""
	align_type = "S" # I guess I'm making my own alignment file now, either P or S based on minimap
	start = 0
	end = 0
	MQ = 0

	def __init__(self, read_name, align_type, start, end, MQ):
		self.read_name = read_name
		self.align_type = align_type
		self.start = start
		self.end = end
		self.MQ = MQ
		

	def __str__(self):
		return "%s\t%s\t%d\t%d\t%d" % (self.read_name, self.align_type, self.start, self.end, self.MQ)


def parseSam(sam_string):
	# Based on sam file
	read_name = sam_string[0]
	flag = int(sam_string[1])
	start = int(sam_string[3])
	MQ = int(sam_string[4])
	cigar = sam_string[5]

	return read_name, flag, start, MQ, cigar

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
		return False
	#####
#####
	


def main():

	# Get sam file
	sam_fh = sys.stdin

	keepSupp = 1
	try:
		keepSupp = int(sys.argv[1])
	except IndexError:
		pass

	# initialize starting variables
	curr_align = None

	# Read through the sam file
	for line in sam_fh:
		read_name, flag, start, MQ, cigar = parseSam(line.split())


		# print("%s\t%d\t%d\t%d" % (read_name, flag, start, MQ))
		align_type = alignType(flag)
		# print(align_type)

		if (align_type):
			# Is an alignment we care about. We are ignoring the reverse alignments for now

			align_length = parseCigar(cigar)
			end = start + align_length

			if (not curr_align):
				# initialize
				curr_align = Alignment(read_name, align_type, start, end, MQ)
			else:
				# Assuming the order of the sam file, if we switch to a new read, the first alignment listed will have an alignment type of "P"
				if (keepSupp):
					if (align_type == "supplementary"):
						# Continue on the current alignment
						if (end > curr_align.end):
							curr_align.end = end

						if (start < curr_align.start):
							curr_align.start = start
					else:
						# Print out the results of the current alignment
						print(curr_align)
						# Start a new alignment
						curr_align = Alignment(read_name, align_type, start, end, MQ)
					#####
				elif (align_type != "supplementary"):
					# Only include the main alignment
					print(curr_align)
					# Start a new alignment
					curr_align = Alignment(read_name, align_type, start, end, MQ)

			#####
		#####
	######
	
	# print the last one
	if (keepSupp or align_type != "supplementary"):
		print(curr_align)
#####


if __name__ == "__main__": main()

