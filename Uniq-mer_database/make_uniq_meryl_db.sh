#!/bin/bash

# This script just makes the databases
# Example 1: ./make_meryl_db.sh chr22_asm 21 unique pos  will produce a database called chr22_asm.sck_pos.meryl with kmer size 21 with unique positions
# Example 2: ./make_meryl_db.sh chr22.sim_reads 21 will produce a database called chr22.sim_reads.meryl with kmer size 21 

# Inputs
prefix=$1		# Database and file names
input_fa=$2		# Fasta to kmerize
k_size=$3		# K-mer size
threads=$4		# Number of threads
memory=$5		# Memory

echo Making a $unique database with prefix $prefix with kmer size $k_size

# Make the unique pos database

main_meryl=$prefix.meryl
echo meryl count k=$k_size $input_fa output $main_meryl
meryl count k=$k_size $input_fa output $main_meryl


# Get unique kmers (value = 1)
equal_to_1_meryl=$prefix.equal_to_1.meryl
meryl equal-to 1 $main_meryl output $equal_to_1_meryl

# Get the positions of the unique forward kmers
to_import_file=$prefix.sck_pos.to_import.txt
meryl-lookup -dump -mers $equal_to_1_meryl -sequence $input_fa -threads $threads -memory $memory | awk -F "\t" '$5==1{print $4"\t"$2}'  > $to_import_file
 # meryl-lookup -dump -mers chr22.asm.sck.meryl -sequence chr22_info/chr22.fasta -threads 8 -memory 20g | awk -F "\t" '$5=1{print $4"\t"$2}' | head

# Make new  postition database
uniq_pos_meryl=$prefix.uniq_pos.meryl
meryl-import -k $k_size -kmers $to_import_file -output $uniq_pos_meryl -threads $threads -memory $memory


# Make kmer list
kmer_list_file=$prefix.kmerlist.txt
meryl-lookup -dump -mers $main_meryl -sequence $input_fa -threads $threads -memory $memory | awk '$5==1{print $4"\t"$5}' > $kmer_list_file





# meryl-import -k 21 -kmers chr22.asm.sck_pos.to_import.txt -output chr22.asm.sck_pos.meryl -threads 8 -memory 20g

echo Done making databases! 




