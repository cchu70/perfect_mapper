#!/bin/bash

# This script is to execute the make meryl database stage of sck analysis script


# Arguments

# inputs
chr=$1
k_size=$2
asm_fa=$3
sim_read_fa=$4
org_read_fa=$5


# Outputs
asm_sck_pos_db=$6
sim_reads_db=$7
org_reads_db=$8
sim_reads_asm_sck=$9
org_reads_asm_sck=$10


# Paths
path_to_scripts=$11



echo Making unique position database 
$path_to_scripts/make_meryl_db.sh $chr.asm $asm_fa $k_size unique position $asm_sck_pos_db
# make unique position datbaase for the assembly

echo Making simulated read database
$path_to_scripts/make_meryl_db.sh $chr.sim_reads $sim_read_fa $k_size $sim_reads_db
# make basic meryl database for the simulated reads

echo Making origin read database
$path_to_scripts/make_meryl_db.sh $chr.org_reads $org_read_fa $k_size $org_reads_db	 
# make basic meryl database for the origin of the sim reads

echo Intersection 
meryl intersect $sim_reads_db $asm_sck_pos_db output $sim_reads_asm_sck	 
# Get the positions and counts of the unique kmers in the reads 

meryl intersect $org_reads_db $asm_sck_pos_db output $org_reads_asm_sck