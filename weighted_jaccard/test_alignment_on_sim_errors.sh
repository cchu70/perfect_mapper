#!/bin/bash

GAGE_A=$1
GAGE_B=$2
which_to_error=$3
which_not_to_error=$4
error_rate=$5
iterations=$6
reads=$7
prefix=$8 # GAGE


# Make new directory and cd into it
error_directory="error_"$error_rate

echo Making $error_directory

cd $error_directory
 echo In directory $error_directory


counter=1

while  [ $counter -le $iterations ] do

	new_fasta_name="$prefix_$which_to_error.err_${error_rate}.v_${iterations}.fasta"

	echo New fasta with error rate $error_rate : $new_fasta_name

	new_split_fasta_name="split.err_${error_rate}_${which_to_error}.fasta"

	echo Write new split fasta : $new_split_fasta_name

	# produce error version of A or B
	if [ $which_to_error = "B" ] 
	then
		which_to_error_fasta=$GAGE_B
		which_not_to_error_fasta=$GAGE_A
	else
		which_to_error_fasta=$GAGE_A
		which_not_to_error_fasta=$GAGE_B
	fi

	echo Introducing error to $which_to_error

	python /data/Phillippy/projects/perfect-polish/scripts/mashmap_postfilter/weighted_jaccard/simulate_error.py $which_to_error_fasta $error_rate $prefix_$which_to_error.err_${error_rate}.v_${iterations}.fasta

	# combine files to create new fasta file
	cat $which_not_to_error_fasta > $new_split_fasta_name
	cat $new_fasta_name >> $new_split_fasta_name

	echo Finished writing new split fasta : $new_split_fasta_name

	# Align the reads onto this new file
	echo Mapping $GAGE_A to $new_split_fasta_name
	sam_A="minimap2.N50_r3k.split.err_${error_rate}_${which_to_error}.aligned_A.sam"
	minimap2 -t12 -a -N50 -r3000 $new_split_fasta_name $GAGE_A -o $sam_A

	echo Wrote bam file to $sam_A

	echo Mapping $GAGE_B to $new_split_fasta_name
	sam_B="minimap2.N50_r3k.split.err_${error_rate}_${which_to_error}.aligned_B.sam"
	minimap2 -t12 -a -N50 -r3000 $new_split_fasta_name $GAGE_B -o $sam_B

	echo Wrote bam file to $sam_B


	# parse through the sam file
	from_A_aligned_A=$(samtools view -F 16 -F 256 -F 2048 $sam_A | awk '$3 == "GAGE_A" {print $0}' | wc -l)
	from_A_aligned_B=$(samtools view -F 16 -F 256 -F 2048 $sam_A | awk '$3 == "GAGE_B" {print $0}' | wc -l)
	from_B_aligned_B=$(samtools view -F 16 -F 256 -F $sam_B | awk '$3 == "GAGE_B" {print $0}' | wc -l)
	from_B_aligned_A=$(samtools view -F 16 -F 256 -F $sam_B | awk '$3 == "GAGE_A" {print $0}' | wc -l)
	echo from_A_aligned_A : $from_A_aligned_A
	echo from_A_aligned_B : $from_A_aligned_B
	echo from_B_aligned_B : $from_B_aligned_B
	echo from_B_aligned_A : $from_B_aligned_A


	# output
	echo ">>>>>>>>>>>>>>>>"
	echo $error_rate"\t"$which_to_error"\t"$from_A_aligned_A"\t"$from_A_aligned_B"\t"$from_B_aligned_B"\t"$from_B_aligned_A"\t"$sam_A"\t"$sam_B
	echo ">>>>>>>>>>>>>>>>"
	((counter++))
done