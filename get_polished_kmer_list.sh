#!/bin/bash

# script to make the uniqe kmer and true kmer list from the POLISHED version of the chrX assembly

module load canu

# make the db
# meryl count k=21 memory=40 threads=12 chrX.fasta output chrX.meryl

# Dump
meryl-lookup -dump -sequence chrX.fasta -mers chrX.meryl -threads 12 -min 1 | awk '$5 >0 {print $4"\t"$5}' > chrX.kmer_list.txt