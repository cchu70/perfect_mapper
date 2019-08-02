# SCript to run through a bunch of sam files and count kmers


import sys



def main():


	sam_fofn = sys.argv[1]
	

	for sam_file in open(sam_fofn, "r"):

		# get the scores for each alignment





if __name__ == "__main__": main()