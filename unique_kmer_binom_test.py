# Script to calculate the P value for the number of unique kmers found on a nanopore read given the error rate and the number of unique kmers

# using /data/chucl/chr22_tests/chr22.org_sim_count_comp.txt

from scipy.stats.distributions import poisson
import sys
import math


def main():
	k = int(sys.argv[1])
	err = float(sys.argv[2])
	org_sim_count_comp = sys.argv[3]
	output_file = sys.argv[4]

	p = (1 - err) ** k
	print(p)

	with open(output_file, "w") as fh:

		for line in open(org_sim_count_comp, "r"):
			org_uniq_cnt = line.split()[1]
			sim_uniq_cnt = line.split()[2]
			if (sim_uniq_cnt != "sim_count"):
				#rv = poisson(org_uniq_cnt)
				#p_val = binom_test(x= int(sim_uniq_cnt), n = int(org_uniq_cnt), p= p)
				#p_val = poisson_probability(int(sim_uniq_cnt), int(org_uniq_cnt), p)
				
				p_val= poisson.pmf(sim_uniq_cnt, int(org_uniq_cnt * p))
				fh.write(line.strip() + "\t" + str(p_val) + "\n")
			else:
				fh.write(line.strip() + "\tp_val\n")
			#####
		#####
#####

# def poisson_probability(actual, mean, p):
#     # naive:   math.exp(-mean) * mean**actual / factorial(actual)

#     # iterative, to keep the components from getting too large or small:
#     p = math.exp(-mean)
#     for i in xrange(actual):
#         p *= mean
#         p /= i+1
#     return p
# #####

if __name__ == "__main__": main()


# Commands for getting the pvalues for sets of reads
#  awk '{if(NR==FNR){keys[$1] = $NF}else{if(keys[$1] && keys[$1] > 0.01){print keys[$1]}}}' ../binom_test/tmp.chr22.unique_kmer_binom_test.txt chr22.k21_s500_none.false.out > chr22.unique_kmer_binom_test.mashmap_false_out.p_val.txt