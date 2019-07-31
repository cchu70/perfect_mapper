# Main script to test the performance of mappers on simulated errors on the reference


import sys
import subprocess


def main():

	# Input files


	GAGE_A = sys.argv[1]
	GAGE_B= sys.argv[2]
	GAGE_A_reads= sys.argv[3]
	GAGE_B_reads= sys.argv[4]

	error_rate_start = float(sys.argv[5])
	error_rate_end = float(sys.argv[6])
	error_rate_step = float(sys.argv[7])

	iterations= int(sys.argv[8])
	
	prefix= sys.argv[9]


	script = sys.argv[10]


	e = error_rate_start
	while e < error_rate_end:

		# Error on A
		pA = subprocess.Popen(["/bin/bash", script, GAGE_A, GAGE_B, 'A', 'B', str(e), str(iterations), GAGE_A_reads, GAGE_B_reads, prefix])

		# Error on B
		pB = subprocess.Popen(["/bin/bash", script, GAGE_A, GAGE_B, 'B', 'A', str(e), str(iterations), GAGE_A_reads, GAGE_B_reads, prefix])

		e += error_rate_step
	#####

	









if __name__ == "__main__": main()