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

# Get the initial scoring columns
python ../../scripts/mashmap_postfilter/unique_kmer_order_score.py $dump_file $map_file $prefix.sort.merge.bed $read_idx $start_idx $end_idx $calc_cores > $scores_out


# Get the mapQ scores
python unique_kmer_mapQ.py $scores_out $mapQ_out $true_reads_bed_file > summary.log
