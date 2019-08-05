#!/bin/bash

# Script outlining how to get the bed file of the true cooridinates of a batch of read names


# Need list of names of interest
read_names=$1
chr_parts=$2 

# Awk script to calculate true cooridinates: this is insanely fast compared to the for loop of grepping
awk -f read_bedfiles.awk chr22.parts.bed chr22.read_names.all.txt > chr22.reads.bed

# Get no name version
cat chr22.prim_true.reads.bed | awk '{print "chr22\t"$2"\t"$3}' > chr22.prim_true.reads.no_name.bed


# Getting the overlap between the bed file and the unique kemr regions
# Looking for LOW overlap

# Look for reads with at MOST some % overlap
bedtools intersect -a chr22.prim_true.reads.no_name.bed -b chr22.sck_positions.merge.bed -f 0.25 -v



# No overlapping regions with unique kmer regions: see where they aligned? and if there were secondary ones
bedtools intersect -a chr22.prim_true.reads.no_name.bed -b chr22.sck_positions.merge.bed -v > chr22.prim_true.reads.NOT_unique.bed


# Prep bedfile to include the original name or information of the gap: put int he 4th column
cat chr22.reads.bed | awk '{print "chr22\t"$2"\t"$3"\t"$1}' > chr22.reads.name_field.bed

# get the fasta file for the original regions
bedtools getfasta -fi chr22.fasta -bed chr22.reads.name_field.bed -name | gzip > chr22.reads.origin.fasta

# Get working verions of canu
export CANU=/data/Phillippy/tools/canu/canu-tip/canu/Linux-amd64/bin/

# See how many unique kmers exist in each of the origin regions
meryl-lookup -existence  -sequence chr22.reads.origin.fasta.gz -mers chr22_unique_pos.meryl -threads 8 -memory 20g > chr22.reads.origin.lookup_exist.txt

# Look at the kmer counts for the simulated reads dump file 

# Plot kmers in sim versus origin read as a dot
# Ideally see just a straight diagonal, but any other discrepancies would fall off the diagonal

# Would see a bunch of dots around the middle ish? like if I plotted a smooth scatter plot of the first histogram I made with the frequency of reads with someunique kmer count

# would be interesting to see the histogram with the origins of the reads

# Just get the kmer count scores
cat chr22.reads.origin.lookup_exist.txt | awk '{print $1"\t"$NF}' > chr22.reads.origin.k_count.txt

# Did the sim reads in parts since the prim false and true exist data already existed
cat week1_tests/test.prim_reads.true.1x.lookup_exist.txt | awk -F '[;\t]' '{print $1"\t"$NF}' > chr22.sim_reads.k_count.txt
cat week1_tests/test.prim_reads.false.1x.lookup_exist.txt | awk -F '[;\t]' '{print $1"\t"$NF}'>> chr22.sim_reads.k_count.txt

# To compare the two: they are not necessarily in order.... (if I just did this directly off the original reads, then I could just go one by one)
# Hash whole set first, then take using smaller subset
awk 'BEGIN{print "group\torigin\tsim"}{if(NR==FNR){k_count[$1]=$2}else {print "all\t"k_count[$1]"\t"$2}}' chr22.reads.origin.k_count.txt chr22.sim_reads.k_count.txt > chr22.k_count.comp_sim_origin.txt

# Plot
Rscript k_count.comp_sim_origin.R chr22.sim_origin_k_count_comp chr22.k_count.comp_sim_origin.txt 10 10

# Do groups? force all the text files to have some group
cat week1_tests/test.prim_reads.true.1x.lookup_exist.txt | awk -F '[;\t]' '{print $1"\t"$NF}' > chr22.sim_reads.prim_true.k_count.txt
cat week1_tests/test.prim_reads.false.1x.lookup_exist.txt | awk -F '[;\t]' '{print $1"\t"$NF}' > chr22.sim_reads.prim_false.k_count.txt

awk 'BEGIN{print "group\torigin\tsim"}{if(NR==FNR){k_count[$1]=$2}else {print "prim_true\t"k_count[$1]"\t"$2}}' chr22.reads.origin.k_count.txt chr22.sim_reads.prim_true.k_count.txt > chr22.k_count.comp_sim_origin.prim_true.txt
awk 'BEGIN{print "group\torigin\tsim"}{if(NR==FNR){k_count[$1]=$2}else {print "prim_false\t"k_count[$1]"\t"$2}}' chr22.reads.origin.k_count.txt chr22.sim_reads.prim_false.k_count.txt > chr22.k_count.comp_sim_origin.prim_false.txt

# The two above lines got the header twice, which messed up the plotter
awk 'BEGIN{print "group\torigin\tsim"}{if(NR==FNR){k_count[$1]=$2}else {print FILENAME"\t"k_count[$1]"\t"$2}}' chr22.reads.origin.k_count.txt chr22.sim_reads.prim_true.k_count.txt chr22.sim_reads.prim_false.k_count.txt > chr22.k_count.comp_sim_origin.prim_groups.txt
# Cat together


# Include the read length as an axis?
# Go back to the exist file
k=21

# Origin
cat week2/chr22.reads.origin.lookup_exist.txt | awk -v k=$k '{print $1"\t"$NF"\t"$2 + k}' > chr22.reads.origin.k_count_len.txt

# Simulated
cat week1_tests/test.prim_reads.true.1x.lookup_exist.txt | awk -F '[;\t]' '{print $1"\t"$NF"\t"$(NF-2)+ 20}' > chr22.sim_reads.k_count_len.txt
cat week1_tests/test.prim_reads.false.1x.lookup_exist.txt | awk -F '[;\t]' '{print $1"\t"$NF"\t"$(NF-2)+ 20}' >> chr22.sim_reads.k_count_len.txt


