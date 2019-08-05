#!/bin/bash

# Script to get the plots for evaluating the loss of uniq-mers and gain of false uniq-mers in simulated nanopore reads


############################################
# Inputs
############################################

# unique meryl DB
uniq_meryl_db=$1

# Simulated nanopore reads
sim_reads_fasta=$2

# Base reference fasta (where the reads were simulated off of)
origin_ref_fasta=$3

# Bedfile with the sim read positions (format)
origin_pos_bedfile=$4

# prefix


############################################
# Script
############################################

# Get the kmer counts on the simulated reads


# Get the kmer counts on the origins of the simulated reads


# Get all the positions of the unique kmers that exist in each read
meryl-lookup -dump -sequence chr22_info/chr22.sim_reads.fasta -mers chr22.asm.sck_pos.meryl -threads 8 -memory 20g | awk '$3=="T"{print $1"\t"$5}' > chr22.sim_reads.asm.sck_pos.dump.txt

# Count the number of true and false uniqmers 
awk 'BEGIN{read="read"; true="true"; false="false"}{if(NR==FNR){start[$4]=$2;end[$4]=$3}else{if ($1 != read){print read"\t"true"\t"false; read = $1; true=0; false=0;} if (start[$1] < $2 && $2 < end[$1]){true=true + 1}else{false=false + 1}}}END{print read"\t"true"\t"false}' chr22.reads.org_pos.bed chr22.sim_reads.asm.sck_pos.dump.txt > chr22.sim_reads.asm.sck_pos.correct_kmer_count.txt



# Count the number fo true and false uniqmers each non-errored version of the simulated reads
bedtools getfasta $origin_pos_bedfile

meryl-lookup -dump -sequence chr22_info/chr22.sim_reads.fasta -mers chr22.asm.sck_pos.meryl -threads 8 -memory 20g | awk '$3=="T"{print $1"\t"$5}' > chr22.origin_reads.asm.sck_pos.dump.txt
cat chr22.origin_reads.asm.sck_pos.dump.txt | awk 'BEGIN { read = ""; count = 0 } { if !read { read = $1 }; if ( $1 != read ) { print read"\t"count; read = $1; count = 1} else { count = count + 1 } }' > origin_reads.true_uniqmer_count.txt



# Compile information

awk ' if ( NR == FNR ){ true_count[$1] = $2 } else { print $0"\t"true_count[$1]} ' origin_reads.true_uniqmer_count.txt chr22.sim_reads.asm.sck_pos.correct_kmer_count.txt > out

############################################
# Plots
############################################

# uniqmer loss plot: the number of the uniqmers (true and false) in the simulated read versus the number that is supposed to be there
cat out | awk '{print $1"\t"($2 + $3)"\t"$4}' > uniqmer_loss.to_plot.txt

Rscript plot

# uniqmer false uniqmer rate
cat out	| awk '{ print $1"\t"( $NF / ( $NF + $(NF - 1) ) ) }' > false_uniqmer_rate.to_plot.txt