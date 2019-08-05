#!/bin/bash

# Get the frequency histogram of the spacing of unique kmers in a genome

# Inputs: Unique kmer meryl database derived from the genome of interest, and the fasta file of the genome of interest

genome=$1 #genome or chromosome name
unique_db=$2
chr_gaps=$3 # bed file of the gaps in the original genome

module load canu
module load bedtools
module load R

# Get the positions of the kmers in the genome
meryl-lookup -dump -sequence $genome.fasta -mers $unique_db -threads 8 -memory 20g | awk '$(NF-4)=="T" {print "chr22\t"$(NF-5)"\t"($(NF-5) + 21)}' > $genome.sck_pos.bed

# Get a merged bedfile of the unique kmer regions
bedtools merge -i $genome.sck_pos.bed > cat $genome.sck_pos.merge.bed

# Get the gap intervals between the unique kmer regions
cat $genome.sck_pos.merge.bed | awk 'BEGIN{left_end=-1;right_front=-1;}{if (left_end < 0) {left_end=($NF-1);} else {right_front=$(NF-1) + 1; print $1"\t"left_end"\t"right_front; left_end=$NF - 1;}}' > $genome.sck_gaps.bed

# Subtract existing gaps in the chromosome
bedtools subtract -a $genome.sck_positions.merge.gaps.bed -b $3 > $genome.sck_gaps.subtract.bed

# Get the sizes of the new intervals (spacing)
cat $genome.sck_gaps.subtract.bed | awk '{print $NF - $(NF - 1)}' > $genome.sck_spacing.txt

# Get frequency histogram
java -jar -Xmx1G /home/rhiea/codes/txtColumnSummary.jar 1 $genome.sck_spacing.txt | awk '{print $1" "$2}' > $genome.sck_spacing.histo

# Make plot
Rscript sck_spacing.plot.R $genome.sck_spacing $genome.sck_spacing.histo 20 10 1000 1000