#!/bin/bash

GAGE_A=$1
GAGE_B=$2
which_to_error=$3
which_not_to_error=$4
error_rate=$5
iterations=$6
GAGE_A_reads=$7
GAGE_B_reads=$8
prefix=$9 # GAGE

out=${prefix}.err_${error_rate}_${which_to_error}.out


module load samtools
module load minimap2

# Make new directory and cd into it
error_directory="error_"$error_rate

mkdir $error_directory
echo ">>>>>>>>>>>>>>>> "Made $error_directory

cd $error_directory
echo ">>>>>>>>>>>>>>>> "In directory $error_directory



if [ ! -f $out ] 
then
	echo Making new file : $out
else
	echo Deleting old file $out
	rm $out
fi


counter=1

while  [ $counter -le $iterations ] 
do

	err_fasta_name="${prefix}.err_${error_rate}_${which_to_error}.v_${counter}.fasta"

	echo ">>>>>>>>>>>>>>>> "New fasta with error rate $error_rate : $new_fasta_name

	new_split_fasta_name="${prefix}_split.err_${error_rate}_${which_to_error}.v_${counter}.fasta"

	echo ">>>>>>>>>>>>>>>> "Write new split fasta : $new_split_fasta_name

	# produce error version of A or B
	if [ $which_to_error = "B" ] 
	then
		which_to_error_fasta=$GAGE_B
		which_not_to_error_fasta=$GAGE_A
	else
		which_to_error_fasta=$GAGE_A
		which_not_to_error_fasta=$GAGE_B
	fi

	echo ">>>>>>>>>>>>>>>> "Introducing error to $which_to_error

	python /data/Phillippy/projects/perfect-polish/scripts/mashmap_postfilter/weighted_jaccard/simulate_error.py $which_to_error_fasta $error_rate $err_fasta_name

	# combine files to create new fasta file
	cat $which_not_to_error_fasta > $new_split_fasta_name
	cat $err_fasta_name >> $new_split_fasta_name

	echo ">>>>>>>>>>>>>>>> "Finished writing new split fasta : $new_split_fasta_name

	# Align the reads onto this new file
	echo ">>>>>>>>>>>>>>>> "Mapping $GAGE_A_reads to $new_split_fasta_name
	sam_A="${prefix}_minimap2.N50_r3k.split.err_${error_rate}_${which_to_error}.v_${counter}.aligned_A.sam"

	minimap2 -t12 -a -N50 -r3000 $new_split_fasta_name $GAGE_A_reads -o $sam_A

	echo ">>>>>>>>>>>>>>>> "Wrote bam file to $sam_A

	echo ">>>>>>>>>>>>>>>> "Mapping $GAGE_B_reads to $new_split_fasta_name
	sam_B="${prefix}_minimap2.N50_r3k.split.err_${error_rate}_${which_to_error}.v_${counter}.aligned_B.sam"
	minimap2 -t12 -a -N50 -r3000 $new_split_fasta_name $GAGE_B_reads -o $sam_B

	echo ">>>>>>>>>>>>>>>> "Wrote bam file to $sam_B


	# parse through the sam file
	from_A_aligned_A=$(samtools view -F 16 -F 256 -F 2048 $sam_A | awk '$3 == "GAGE_A" {print $0}' | wc -l)
	from_A_aligned_B=$(samtools view -F 16 -F 256 -F 2048 $sam_A | awk '$3 == "GAGE_B" {print $0}' | wc -l)
	from_B_aligned_B=$(samtools view -F 16 -F 256 -F 2048 $sam_B | awk '$3 == "GAGE_B" {print $0}' | wc -l)
	from_B_aligned_A=$(samtools view -F 16 -F 256 -F 2048 $sam_B | awk '$3 == "GAGE_A" {print $0}' | wc -l)

	echo ">>>>>>>>>>>>>>>> "from_A_aligned_A : $from_A_aligned_A
	echo ">>>>>>>>>>>>>>>> "from_A_aligned_B : $from_A_aligned_B
	echo ">>>>>>>>>>>>>>>> "from_B_aligned_B : $from_B_aligned_B
	echo ">>>>>>>>>>>>>>>> "from_B_aligned_A : $from_B_aligned_A


	# output
	echo -e $error_rate"\t"$which_to_error"\t"$from_A_aligned_A"\t"$from_A_aligned_B"\t"$from_B_aligned_B"\t"$from_B_aligned_A"\t"$sam_A"\t"$sam_B >> $out
	((counter++))
done

echo Finished




