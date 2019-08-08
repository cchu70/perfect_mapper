#!/bin/bash

# Script to retrieve the new ground truth using the best alignmnet of the non-simulated reads onto the unpolished assembly


# Inputs


# # fw.list is the list of the simulated reads
# cat fw.list | awk -F '[_\t]' '{print "chrX_fixedBionanoSV_centromereV4_racon_patch139_arrow_arrow\t"$2"\t"$3"\t"$0}' > fw.chrX.sim_reads.origin.bed

# # retrieve the original sequence of these reads
# bedtools getfasta -fi chrX.fasta -bed fw.chrX.sim_reads.origin.bed -name -fo chrX.orig_reads.fasta


# Align these origin reads onto the prepolished chrX
minimap2 -t12 -a -N50 -r3000 chrX.prepolished.fasta chrX.orig_reads.fasta -o chrX.prepolished.orig_reads.minimap2.N50_r3k.sam

# filter the samtools alignment to rempve split and only consider primary reads. This is so each read gets only one alignment. Return the start position of the read
samtools -F 2048 -F 256 -F 4 -F 16 view chrX.prepolished.orig_reads.minimap2.N50_r3k.sam | awk '$2 == 0 {print $1"\t"$4}' > chrX.prepolished.orig_reads.ground_truth_pos.txt
 d
