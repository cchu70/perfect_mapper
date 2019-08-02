#!/bin/bash

# Script to run weighted_jaccard_count_plain_sam_input.py to test performance on GAGE locus


prefix=$1
err_start=$2
err_end=$3
err_step=$4

v_start=$5
v_end=$6

i=$err_start
v=$v_start

while [ $i -le $err_end ] 
do
	while [ $v -le $v_end ]
	do
		echo "python ../../scripts/mashmap_postfilter/weighted_jaccard/weighted_jaccard_count_plain_sam_input.py GAGE_A.sim_reads.fasta error_${i}/${prefix}.err_${i}_A.v_1.fasta error_${i}/${prefix}_minimap2.N50_r3k.split.err_${i}_A.v_1.aligned_A.sam GAGE.kmerlist.txt 21 GAGE_A GAGE_A ${i}"
		((v++))
	done
	((i = i + err_step))
done