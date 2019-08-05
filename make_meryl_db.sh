#1/bin/bash

# This script just makes the databases
# Example 1: ./make_meryl_db.sh chr22_asm 21 unique pos  will produce a database called chr22_asm.sck_pos.meryl with kmer size 21 with unique positions
# Example 2: ./make_meryl_db.sh chr22.sim_reads 21 will produce a database called chr22.sim_reads.meryl with kmer size 21 

# Inputs
prefix=$1
input_fa=$2
output_db=$3
k_size=$3
unique=$4 # set to 1 if yes, 0 for no
pos=$5
threads=$6
memory=$7

echo Making a $unique database with prefix $prefix with kmer size $k_size

# Make the unique pos database
echo meryl count k=$k_size $input_fa output $prefix.meryl 
meryl count k=$k_size $input_fa output $prefix.meryl

if [ $unique == "unique" ] # Get unique kmers only
then
	# Get unique kmers (value = 1)
	meryl equal-to 1 $prefix.meryl output $output
else if [$pos == "position"]
then
	# Get the positions of the unique forward kmers
	meryl-lookup -dump -mers $prefix.sck.meryl -sequence $input_fa -threads $threads -memory $memory > $prefix.sck.lookup_dump.txt | awk -F "\t" '$5=1{print $4"\t"$2}'  > $prefix.sck_pos.to_import.txt
	 # meryl-lookup -dump -mers chr22.asm.sck.meryl -sequence chr22_info/chr22.fasta -threads 8 -memory 20g | awk -F "\t" '$5=1{print $4"\t"$2}' | head

	# Make new database
	meryl-import -k $k_size -kmers $prefix.sck_pos.to_import.txt -output $output_db -threads $threads -memory $memory
	# meryl-import -k 21 -kmers chr22.asm.sck_pos.to_import.txt -output chr22.asm.sck_pos.meryl -threads 8 -memory 20g
fi

echo Done making database! 

# # Find the unique kmers in the reads from the unique database

# # Get fasta files
# sim_reads.fa=$something

# ./get_fasta_with_bed readNames chr_fa_file 

# meryl-lookup -existence -sequence sim_reads.fa -mers ${chr}.asm.sck_pos.meryl> ${chr}.sim_reads.sck_count.txt
# meryl-lookup -existence -sequence org_reads.fa -mers ${chr}.asm.sck_pos.meryl> ${chr}.org_reads.sck_count.txt



