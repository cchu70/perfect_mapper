#!/bin/bash

# Script to run weighted_jaccard_count_plain_sam_input.py to test performance on GAGE locus


prefix=$1
err_start=$2
err_end=$3
err_step=$4

v_start=$5
v_end=$6

which_part=$7

i=$err_start
v=$v_start

while [ $i -le $err_end ] 
do
	while [ $v -le $v_end ]
	do
		python ../../scripts/mashmap_postfilter/weighted_jaccard/weighted_jaccard_count_plain_sam_input.py GAGE_${which_part}.sim_reads.fasta error_${i}/${prefix}.err_${i}_${which_part}.v_1.fasta error_${i}/${prefix}_minimap2.N50_r3k.split.err_${i}_${which_part}.v_1.aligned_${which_part}.sam GAGE.kmerlist.txt 21 GAGE_${which_part} GAGE_${which_part} ${i}
		((v++))
	done
	((i = i + err_step))
done