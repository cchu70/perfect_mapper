#!bin/bash


read_fasta="/data/Phillippy/projects/perfect-polish/chm13_chrX/chrX-a02-s10.simulated.fw.fasta"
ref_fasta="/data/Phillippy/projects/perfect-polish/chm13_chrX/chrX.renamed.fasta"
align_file="chrX.unpolished.minimap2.default.bam"



module load samtools
module load canu


# Samtools get only the forward alignments

samtools view -F 2048 -F 16 -F 4

# Get the true positions of the reads

awk -f script that got the origin stuff > chrX-a02-s10.org_pos.bed


# Get the ground truth on the bam file
samtools view $align_file | python /data/Phillippy/projects/perfect-polish/scripts/mashmap_postfilter/get_multiple_alignments_from_bam.py chrX-a02-s10.org_pos.bed > chrX.unpolished.minimap2.default.mult_align.ground_truth.sam

# Make the unique kmer list

meryl count $ref_fasta output chrX.unpolished.meryl

meryl-lookup -dump -mers chrX.unpolished.meryl -seq $ref_fasta | awk '{print only the forward kmers' > chrX.unique_kmer_list.txt


# run weighted jaccard

python /data/Phillippy/projects/perfect-polish/scripts/mashmap_postfilter/weighted_jaccard/weighted_jaccard_main.py $read_fasta $ref_fasta chrX.unpolished.minimap2.default.mult_align.ground_truth.sam chrX.unique_kmer_list.txt 1 16 21
