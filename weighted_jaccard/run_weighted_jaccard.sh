#!/bin/bash


# Example script to evaluate what scheme weights for uniq-mers produces the best alignments


# Inputs

sim_reads_fa=$1				# Fasta file with the simulated reads used to produce mapping file
ref_fa=$2					# fasta file with single header of the chromosome or ref the reads are mapped to
map_file=$3				# Alignment or mapping file containing the following information: Read name, read start index, target start index, target end index. depending on type, will also have Primary or Secondary flags (sam), or % idy (mashmap) (Ex. chrX.prepolished.chrX-a02-s10.simulated.mashmap.out)
map_file_type=$4			# Ex. mashmap [or sam]
ground_truth_bedfile=$5 	# Bedfile in the form "read_name<tab>real_start_index<tab>real_end_index" for each of the simulated reads mapped
uniqmer_list=$6				# List of the uniq-mers from the target
k_size=$7 					# Size of the kmer used
scheme_start=$8				# Lowest weight to test on uniq-mer weighted Jaccard scoring
scheme_end=$9				# highest (inclusive) weight to test on uniq-mer weighted jaccard scoring
scheme_step=${10}				# Step sizes of schemes to test between the start and end 
prefix=${11}					# prefix for the output files



# Filtering and reformatting the alignments to only get reads with multiple alignments and relevant information. Additionally appending the ground truth at the end of the line

if [ $map_file_type == 'mashmap' ]: then
	python path/to/get_multiple_alignments_from_mashmap.py $map_file $ground_truth_bedfile > $prefix.mult_align.ground_truth.out

else:
	python path/to/get_multiple_alignments_from_bam.py $map_file $ground_truth_bedfile > $prefix.mult_align.ground_truth.out
fi

# Count the shared and non-shared uniq-mers and non-uniq-mers in the read and the aligned regions
python path/to/weighted_jaccard/weighted_jaccard_count.py $sim_reads_fa $ref_fa ${prefix}.mult_align.ground_truth.outmashmap $uniqmer_list 21 > $prefix.mult_align.ground_truth.k_counts.txt

# Score using different weights for the uniq-mers using the k-mer counts above
python path/to/weighted_jaccard/weighted_jaccard_scheme_score.py $prefix.mult_align.ground_truth.k_counts.txt $scheme_start $scheme_end $scheme_step > $prefix.mult_align.ground_truth.k_counts.sch_scores.${scheme_start}_${scheme_end}_${scheme_step}.txt

# Test the performance
python path/to/weighted_jaccard/weighted_jaccard_eval_main.py $prefix.mult_align.ground_truth.k_counts.sch_scores.${scheme_start}_${scheme_end}_${scheme_step}.txt $prefix.mult_align.ground_truth.k_counts.sch_scores.${scheme_start}_${scheme_end}_${scheme_step}