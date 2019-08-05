#!/bin/bash

# This is the main script to run the full analysis of 
# unique single copy kmers (or sck's) in nanopore reads, 
# chromosome assemblies, and its performance comapred to 
# current alignment strategies


chr=$1 # Input the name of the chromosome to be used. Note that this will be used to reference existing files in your current working directory need (listed below)
k_size=$2
threads=$3
mem=$4
path_to_initial_files=$5
path_to_scripts=$6
bam=$7 # Alignment file from mapping the reads (read_fa) to the assembly (asm_fa)

echo $chr, $k_size, $threads, $mem

# Script for initializing
# Load modules
# set paths, variables 

#############################################
# Initialize names
#############################################

# Input files
asm_fa=$path_to_initial_files/$chr.fasta
asm_part_bed=$path_to_initial_files/$chr.parts.bed
asm_gaps_bed=$path_to_initial_files/$chr.gaps.bed


# Load Modules
module load bedtools
module load samtools
module load R
module load canu
#export CANU=/data/Phillippy/tools/canu/canu-tip/canu/Linux-amd64/bin/


#############################################
# Make fasta files
echo Making fasta files
#############################################

# File names
# Inputs
sim_read_names=${chr}.sim_reads.names.txt
org_read_bed=$chr.reads.org_pos.bed

# Outputs
sim_read_fa=$path_to_initial_files/$chr.sim_reads.fasta # might be simulated later
org_read_fa=$chr.org_reads.fasta

sbatch --cpus-per-task $threads --memory $mem --wait  $path_to_scripts/make_fastas_main.sh $path_to_scripts $sim_read_names $org_read_bed $sim_read_fa $org_read_fa 


#############################################
# Make the meryl databases
echo Making Meryl Databases
#############################################

# File names outputs 
asm_sck_pos_db=$chr.asm.sck_pos.meryl
sim_reads_db=$chr.sim_reads.meryl
org_reads_db=$chr.org_reads.meryl
sim_reads_asm_sck=$chr.sim_reads.asm_sck.meryl
org_reads_asm_sck=$chr.org_reads.asm_sck.meryl

# # inputs
# chr=$1
# k_size=$2
# asm_fa=$3
# sim_read_fa=$4
# org_read_fa=$5


# # Outputs
# asm_sck_pos_db=$6
# sim_reads_db=$7
# org_reads_db=$8
# sim_reads_asm_sck=$9
# org_reads_asm_sck=$10

sbatch --cpus-per-task $threads --memory $mem --wait $path_to_scripts/make_meryl_db_main.sh $chr $k_size $asm_fa $sim_read_fa $org_read_fa $asm_sck_pos_db $sim_reads_db $org_reads_db $sim_reads_asm_sck $org_reads_asm_sck


#############################################
# Spacing of unique kmer regions
#############################################
# Get bedfile
	# get_batch_reads._bedfile.sh
	# sck_spacing.awk
	# sck_spacing.sh
# merge regions
# get the gaps
# subtract existing gaps
# Histogram it

sbatch --cpus-per-task $threads --memory $mem ./sck_spacing_main.sh $chr.asm.sck_pos.meryl 



#############################################
# For all reads histograms: the frequency of unique kmers on each read
#############################################

# exist read kmer in database -> unique kmer count per read
# Frequency histogram
# sck_dist.sh

sbatch --cpus-per-task $threads --memory ./sck_histograms_main.sh 


#############################################
# Split Bam
#############################################
# Split all the bam into prim and all secondary
	# get_bam_from_read_id.sh
	# minimap2_prim_sec_read_compare.sh
	# prim_false.sec_true.awk
	# SAM2FA>awk
	# split_read_by_start_index.sh
# evaluate true or false
# histograms all of them like above
# for select few, make kmer dot plots of the true and alignment regions 
	# single_read.k_dotplot.sh and single_read.k_dotplot.main.sh

sbatch --cpus-per-task $threads --memory ./sck_split_main.sh 

#############################################
# Read length vs. unique kmer count dot plot
#############################################
# readLen_uniquekmer.sh
sbatch --cpus-per-task $threads --memory ./readLen_sck_main.sh 




#############################################
# Read simulated versus true kmer count
#############################################
# get_batch_reads_bedfile.sh
# some Rplots

sbatch --cpus-per-task $threads --memory ./read_org_sim_comp_main.sh 



# Plotting files
# readLen_sck_comp.plot.R

# All the files
	# cal_k_coor.awk								This script parses the file containing the read name,	location of kmer in the	read, and location of kmer in the referenc
	# compare_prim_sec.sh 							Create the histograms for both the primary and secondary reads
	# Gepard-1.40.jar								Skip
	# get_k_coor.sh 								This is just a script to submit to sbatch to get the coordinates of kmers
	# k_coor.dot_plot.R 							This script is to plot unique kmer positions compared with the reference and nanopore reads
	# k_count.comp_sim_origin.R 					Script to plot the differences between the number of unique kmers in the simulated and the origin reads
	# k_count.comp_sim_origin.smooth_scatter.R 		Script to plot the differences between the number of unique kmers in the simulated and the origin reads; Attempting a smooth scatter plot
	# k_count_len.comp_sim_origin.R 				Script to plot the differences between the number of unique kmers in the simulated and the origin reads; 3D plot including the read length
	# make_hist.prim_reads.true.R 					mmm don't know, but I need to come up with an all purpose plotting script anyways
	# make_hist.sh 									This script takes in a lookup table from meryl lookup and converts it into a histogram with frequencies corresponding to the number reads with a certain number of (unique) kmers
	# minimap2_ftest.sh 							Testing if changing the f value will help minimap alignment
	# overlay_hist.prim_reads.true.R 				R script to plot the frequency of unique kmers on reads that are separated into groups
	# overlay_hist.R 								I think I just linked this or copied from Arang
	# plot_k_coor.sh 								This script is to plot the dot plot of unique kmers between the reads and the reference; pretty useless right now
	# prim_false.sec_true.awk 						Purpose is to split the primary false reads into sets where they have a corresponding secondary read that was marked as true
	# prim_false.sec_true.sh 						Script to sort out false prim reads based on if they exist in the secondary reads that are true
	# read_bedfiles.awk 							SCript to create bedfile of the true start and end positions of a read; read names in the format "chr00_part00_<part start (0-based)>_<part end (1 base)>
	# readLen_sck.comp.plot.R 						Plot dot plot of readLen vs unique kmer counts, colorized based on which sequence came from which group
	# readLen_sck_comp.sh 							For each read in the output of meryl-lookup -existence, output a tab delimited file with Col1: group name; col2: read length; col3: unique kmers
	# sck_dist.sh 									This script produces the kmer histogram of the intersection between the single copy kmers in the reference genome to the kmer counts in the raw reads to gauge the distribution of unique kmers throughout the genome
	# sck_spacing.dot_plot.R 						Plot histogram of the spacing between unique kmer regions
	# sck_spacing.plot.R 							Plot histogram of the spacing between unique kmer regions (line plot)
	# sck_spacing.sh 								Get the frequency histogram of the spacing of unique kmers in a genome; Inputs: Unique kmer meryl database derived from the genome of interest, and the fasta file of the genome of interest
	# single_read.k_dotplot.main.sh 				Script to run multiple reads to view unique kmer dotplots
	# single_read.k_dotplot.sh 						Script to make a dot plot for a single read with multiple alignments
	# split_reads_correct.1xShiftLeft.awk 			Splits reads into separate fasta files based on if the aligned region is in a tolerable range (+/- 50%)
	# split_reads_correct.2x.awk 					Splits reads into separate fasta files based on if the aligned region is in a tolerable range (+/- 99.99999%)


# Real files I need to make
	# This one
	# A general no-name Rscript for dot plots, line plots, histograms, or one for each (like R.dot_plot) and can set own parameters, or just if statements to select type fo graph
	# Generalized awk scripts for splitting: variable names and outputs
	# Separate out the smaller awk scripts into individual files with good comments and a descriptive name
	# Scripts for each stage
	# a submit batch script 
	# Some file structure (like for each stage, they get their own folder, handling images, etc.)









