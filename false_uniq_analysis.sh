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
sim_read_dump_file=${prefix}.sim_reads.uniqmers.txt
meryl-lookup -dump -sequence $sim_reads_fasta -mers $uniq_meryl_db -threads 8 -memory 20g | awk '$5 > 0 {print $1"\t"$5}' > $sim_read_dump_file

# Count the number of true and false uniqmers 
sim_read_true_false_uniqmer_count=${prefix}.sim_reads.sim_read_true_false_uniqmer_count.txt
awk 'BEGIN{read="read"; true="true"; false="false"}{if(NR==FNR){start[$4]=$2;end[$4]=$3}else{if ($1 != read){print read"\t"true"\t"false; read = $1; true=0; false=0;} if (start[$1] < $2 && $2 < end[$1]){true=true + 1}else{false=false + 1}}}END{print read"\t"true"\t"false}' $origin_pos_bedfile $sim_read_dump_file > $sim_read_true_false_uniqmer_count



# Count the number fo true and false uniqmers each non-errored version of the simulated reads
origin_reads_fasta=${prefix}.origin_reads.fasta
bedtools getfasta $origin_pos_bedfile -name > $origin_reads_fasta


# Get the kmers in the origin reads
origin_read_dump_file=${prefix}.origin_reads.uniqmers.txt
meryl-lookup -dump -sequence chr22_info/chr22.sim_reads.fasta -mers chr22.asm.sck_pos.meryl -threads 8 -memory 20g | awk '$5 > 0{print $1"\t"$5}' > $origin_read_dump_file

origin_read_true_uniqmer_count=$${prefix}.sim_reads.sim_read_true_uniqmer_count.txt
cat $origin_read_dump_file | awk 'BEGIN { read = ""; count = 0 } { if (read) { if ( $1 != read ) { print read"\t"count; read = $1; count = 1} else { count = count + 1 } } else { read = $1 } }' > $origin_read_true_uniqmer_count

# Compile information
reads_compiled_uniqmer_counts_outfile=${prefix}.reads_compiled_uniqmer_counts.txt
awk ' if ( NR == FNR ){ true_count[$1] = $2 } else { print $0"\t"true_count[$1]} ' $origin_read_true_uniqmer_count $sim_read_true_false_uniqmer_count > $reads_compiled_uniqmer_counts_outfile

############################################
# Plots
############################################

# uniqmer loss plot: the number of the uniqmers (true and false) in the simulated read versus the number that is supposed to be there
cat out | awk '{print $1"\t"($2 + $3)"\t"$4}' > uniqmer_loss.to_plot.txt

Rscript plot

# uniqmer false uniqmer rate
cat out	| awk '{ print $1"\t"( $NF / ( $NF + $(NF - 1) ) ) }' > false_uniqmer_rate.to_plot.txt