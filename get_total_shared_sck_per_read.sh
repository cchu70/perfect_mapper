#!/bin/bash
# Script to get the merged indices of each read's alignments on the reference

paf_file=$1
prefix=$2

module load bedtools


cat $paf_file | awk '{print $1"\t"$8"\t"$9}' > $prefix.bed


# Sort the start and end for each read
bedtools sort -i $prefix.bed > $prefix.sort.bed

# merge
bedtools merge -i $prefix.sort.bed > $prefix.sort.merge.bed
