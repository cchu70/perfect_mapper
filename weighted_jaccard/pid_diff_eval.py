# Script to take the output of the calc_percent_identity.py file and filter and calculate the difference between a primary and secondary alignment's pids, and if they are correct or not


import sys

def main():
	# Get the pid file with the alignment information
	pid_file = sys.argv[1]

	# track current read

	curr_read_name = None
	curr_prim_pid = 0
	curr_sec_pid = 0

	curr_prim_is_true = False
	got_second = False

	for line in open(pid_file, "r"):

		# Get if it is primary 
		# get it's ground truth
		# get it's second best alignment OR if primary is false and this one is false, keep looking until you get to the true alignment, otherwise skip the whole read
		# calculate the difference
		# If minimap is wrong, put it into the incorrect file, else in the correct file. 

		read_name, map_truth, ref_start, ref_end, ground_truth, pid = line.split()

		if curr_read_name:

			if (curr_read_name != read_name):
				# analyze the current read

				if got_second:
					# there is a second one to compare to such that between the prim and sec, at least one of them have a ground truth value of True
					# Subtract the difference in pid between the True and False alignment (not necessary the primary and secondary)
					if curr_prim_is_true:
						diff = curr_prim_pid - curr_sec_pid
						print("correct\t%d" % diff)
					else:
						diff = curr_sec_pid - curr_prim_pid
						print("incorrect\t%d" % diff)
				#####

				# else continue


				# re-init
				curr_read_name = read_name	
				curr_prim_pid = 0
				curr_sec_pid = 0

				curr_prim_is_true = False
				got_second = False
			#####
		else:
			# initialize read name
			curr_read_name = read_name
		#####

		if map_truth == "P":
			curr_prim_pid = pid
			if(ground_truth == "True"):
				# Just take the second alignment
				curr_prim_is_true = True
		else:
			# check if the primary one is already true. If so, just take this one and skip the rest

			# else, keep looking unil you find a true alignment

			if not got_second:
				if(curr_prim_is_true):
					curr_sec_pid = float(pid)
					got_second = True
				else:
					# prim is not true. check if this one is true
					if ground_truth == "True":
						curr_sec_pid = float(pid)
						got_second = True
					#####
				#####
			#####
		#####




if __name__ == "__main__": main()