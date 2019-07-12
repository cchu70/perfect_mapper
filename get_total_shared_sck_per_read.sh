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


# ../../scripts/_submit_norm.sh 20 20g minimap_supp_merge_sck_count get_minimap_paf_sck_counts.sh "../../chr22_info/chr22.asm.sck_pos.intersect.sim_reads.dump.txt ../../chr22_info/chr22.minimap2_N50_30kb.supp_merged.sam chr22.minimap2_N50_30kb.supp_merged.ref_align.sort.bed 0 2 3 count_order_plus1_only,count_shared_sck chr22.minimap2_N50_30kb.supp_merged.total_sck_per_read.order.sck_counts.txt"