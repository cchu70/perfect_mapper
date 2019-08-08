#!/bin/bash

# Script to make a dot plot for a single read with multiple alignments

read_name=$1
chr=$2
sam=$3
unique_pos_db=$4
k=$5
read_bed=$6
reads_fa=$7


# test if the read originates froma region with unique kmers

# grep bed file?
# samtools view $sam | awk '{print $1}' | sort -u


# extract the single read fasta

# echo $read_name | awk -F "_" '{print $1"\t"$3"\t"$4}' > $read_name.bed
#WRONG ^^

grep $read_name $read_bed | awk -v chr=$chr '{print chr"\t"$2"\t"$3}' > $read_name.bed


# Simulated read
cat $reads_fa | awk -v read_name=$read_name '$0 ~ read_name{p++;print;next} /^>/{p=0} p' > ${read_name}.fasta

# REAL REGION
bedtools getfasta -fi $chr.fasta -bed $read_name.bed > $read_name.real.fasta

read_len=$(cat ${read_name}.bed | awk '{print $NF - $(NF-1)}')

echo single read k dotplot for : $read_name >> $read_name.log

echo Read length: $read_len >> $read_name.log


# UNRELIABLE for split reads
# grep $read_name $sam | awk -v read_name=$read_name 'BEGIN{print ">"read_name;}{if ($10 != "*"){print $10}}' > ${read_name}.fasta

meryl-lookup -dump -sequence ${read_name}.fasta -mers $unique_pos_db -threads 8 -memory 20g | awk '($5 + $7) > 0 {print $0}' > ${read_name}.unique_pos.dump.txt

# Originally was assuming forward mers only... not sure how to handle this
cat ${read_name}.unique_pos.dump.txt| awk '{if ($5 > 0) {print $4, $2}else {print $6, $2}}' | awk '{print $1" "$2 + 1}' > ${read_name}.unique_pos_plus1.to_import.txt

echo Unique kmer count: >>  $read_name.log
wc -l ${read_name}.unique_pos_plus1.to_import.txt >> $read_name.log


# this only includes the kmers in the read that was also in the database, which are all unique kmers


# Make the single read database
meryl-import -kmers ${read_name}.unique_pos_plus1.to_import.txt -k $k -output  ${read_name}.unique_pos_plus1.meryl -threads 8 -memory 20g



grep $read_name $sam | awk -v chr=$chr -v read_len=$read_len '{print chr"\t"($4 - 50)"\t"($4 + read_len + 50)}' > ${read_name}.ref_regions.bed

bedtools merge -i ${read_name}.ref_regions.bed > ${read_name}.ref_regions.merge.bed # maybe not

##### Getting the reduced fasta file 
# reduced fasta, including the original
cat $read_name.real.fasta > ${read_name}.ref_regions.fasta
bedtools getfasta -fi $chr.fasta -bed ${read_name}.ref_regions.merge.bed >> ${read_name}.ref_regions.fasta


 # getting the sequence

# Extract all the alignment regions the read was mapped to from the ref or assembly

# # get combined fasta file of these regions with "N"'s to separate them 
# cat ${read_name}.ref_regions.fasta | awk '!/>/{print $0}' | awk 'ORS="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"' >> ${read_name}.ref_regions.combine.fasta
# # NOPE ^^





meryl-lookup -dump -sequence ${read_name}.ref_regions.fasta -mers  ${read_name}.unique_pos_plus1.meryl/ -threads 8 -memory 20g >  ${read_name}.unique_pos_plus1.ref_regions.dump.txt 
# Make coordinates
 cat  ${read_name}.unique_pos_plus1.ref_regions.dump.txt | awk '$5>0{print $1"\t"$2"\t"$5}' >  ${read_name}.unique_pos_plus1.ref_regions.k_coor.TMP

#  cat  ${read_name}.unique_pos_plus1.ref_regions.k_coor.TMP | awk -F '[:\t-]' '{print $2+$4"\t"$5-1}' >  ${read_name}.unique_pos_plus1.ref_regions.k_coor.txt

# # Above method, since I am actually getting hte genomic regions, the grap is just look like vertical lines


# # To do the combined version...
# meryl-lookup -dump -sequence  ${read_name}.ref_regions.combine.fasta -mers  ${read_name}.unique_pos_plus1.meryl/ >  ${read_name}.unique_pos_plus1.ref_regions.combine.dump.txt -threads 8 -memory 20g

#  cat  ${read_name}.unique_pos_plus1.ref_regions.combine.dump.txt | awk 'BEGIN{print "ref\tread"}{if ($5>0){print $2"\t"$5}}' >  ${read_name}.unique_pos_plus1.ref_regions.combine.k_coor.txt

# Rscript  k_coor.dot_plot.R  ${read_name}.unique_pos_plus1.ref_regions.combine  ${read_name}.unique_pos_plus1.ref_regions.combine.k_coor.txt 10 20



# # Or just extract from an earlier dump? from 
# # meryl-lookup -dump -sequence chr22-a02-s10.reads.fasta -mers chr22_unique_pos.meryl/ -threads 8 -memory 20g > chr22.reads.nuqiue_pos.dump.txt
# # but this file was gigantic (over 40 G), so not feasible

# # Trying color separation? because I do not know which is which

# cat  ${read_name}.unique_pos_plus1.ref_regions.k_coor.TMP | awk 'BEGIN{print "group\tref\tread"}{print $1"\t"($2)"\t"$3-1}' >  ${read_name}.unique_pos_plus1.ref_regions.k_coor.groups.txt 

# Rscript k_coor.dot_plot.R  ${read_name}.unique_pos_plus1.ref_regions.k_coor.groups  ${read_name}.unique_pos_plus1.ref_regions.k_coor.groups.txt 20 10


# Or find some funky way to split up the parts...

cat  ${read_name}.unique_pos_plus1.ref_regions.k_coor.TMP | awk -v read_name=$read_name 'BEGIN{print "group\tref\"read_name; start = 0; count = 0; curr=""}{if ($1==curr){count = $2;} else {curr = $1; start = start + count;} print $1"\t"($2 + start)"\t"$3-1;}' >  ${read_name}.unique_pos_plus1.ref_regions.k_coor.sketchy_group.txt
Rscript k_coor.dot_plot.R  ${read_name}.unique_pos_plus1.ref_regions.k_coor.sketchy_groups  ${read_name}.unique_pos_plus1.ref_regions.k_coor.sketchy_group.txt 10 10




 # Wow that actually worked...

