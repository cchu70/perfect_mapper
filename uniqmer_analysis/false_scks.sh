# Scripts for finding false unique kmers


# comparing the total number of unique kmers




# Get all the positions of the unique kmers that exist in each read
meryl-lookup -dump -sequence chr22_info/chr22.sim_reads.fasta -mers chr22.asm.sck_pos.meryl -threads 8 -memory 20g | awk '$3=="T"{print $0}' > chr22.sim_reads.asm.sck_pos.dump.txt

# Probably include an awk script at the end to further parse? OR you just use the plain coordinate version of the read fasta file? Nah, that is only for the origin

# have the bed file chr22.reads.org_pos.bed: if the kmer is within range
awk '{if(NR==FNR){start[$4]=$2;end[$4]=$3}else{if (start[$1] < $5 && $5 < end[$1]){print $0}}}' chr22.reads.org_pos.bed chr22.sim_reads.asm.sck_pos.dump.txt 

# for the sake of space, just keep the read name and the position of the kmer in the assembly
meryl-lookup -dump -sequence chr22_info/chr22.sim_reads.fasta -mers chr22.asm.sck_pos.meryl -threads 8 -memory 20g | awk '$3=="T"{print $1"\t"$5}' > chr22.sim_reads.asm.sck_pos.dump.txt
awk '{if(NR==FNR){start[$4]=$2;end[$4]=$3}else{if (start[$1] < $2 && $2 < end[$1]){print $0}}}' chr22.reads.org_pos.bed chr22.sim_reads.asm.sck_pos.dump.txt 

# file with the original kmer counts
# chr22.k_count.comp_sim_origin.txt

# group	sim	origin
# all	2878	5180
# all	5626	9459
# all	485	998
# all	2633	5355
# all	1540	2883
# all	3059	5510
# all	4092	7446
# all	4807	8867
# all	3988	6918

# first sim read kcount chr22.sim_reads.k_count.txt
# chr22_part01_8_9885	2878
# chr22_part01_325_19686	5626
# chr22_part01_434_2262	485
# chr22_part01_542_10612	2633
# chr22_part01_712_5737	1540
# chr22_part01_1649_12304	3059
# chr22_part01_1693_16774	4092
# chr22_part01_2215_20824	4807
# chr22_part01_2829_16998	3988
# chr22_part01_3003_8864	1554

# So just modify the awk statement to count the number of unique kmers within region

awk 'BEGIN{read="read"; true="true"; false="false"}{if(NR==FNR){start[$4]=$2;end[$4]=$3}else{if ($1 != read){print read"\t"true"\t"false; read = $1; true=0; false=0;} if (start[$1] < $2 && $2 < end[$1]){true=true + 1}else{false=false + 1}}}END{print read"\t"true"\t"false}' chr22.reads.org_pos.bed chr22.sim_reads.asm.sck_pos.dump.txt > chr22.sim_reads.asm.sck_pos.correct_kmer_count.txt


# Compare this text with the sim_reads.k_count.txt

awk 'BEGIN{print "read\ttrue\tactual"}{if(NR==FNR){true_count[$1]=$2;}else{print $1"\t"true_count[$1]"\t"$2}}' chr22.sim_reads.asm.sck_pos.correct_kmer_count.txt week2/chr22.sim_reads.k_count.txt

# May need to check that I am only looking at forward strands
meryl-lookup -existence -sequence chr22_info/chr22.sim_reads.fasta -mers chr22.asm.sck_pos.meryl -threads 8 -memory 20g | awk '{print $1"\t"$NF}' > chr22.sim_reads.asm.sck_pos.meryl.exist.txt
meryl-lookup -existence -sequence chr22.org_reads.fasta -mers chr22.asm.sck_pos.meryl -threads 8 -memory 20g | awk '{print $1"\t"$NF}' > chr22.org_reads.asm.sck_pos.meryl.exist.txt
 # Now the numbers look right? maybe because before it was based off split prim true and false
 # Compare again

awk 'BEGIN{print "read\torg_count\tsim_count\tsim_true\tsim_false"}{if(NR==FNR){org_count[$1]=$2;}else{ if ($1 ~ /chr/){print $1, org_count[$1]"\t"($2 + $3)"\t"$2"\t"$3}}}' chr22.org_reads.asm.sck_pos.meryl.exist.rename.txt chr22.sim_reads.asm.sck_pos.correct_kmer_count.txt > chr22.sck_count_comp.txt

# Remaking dotplot
cat chr22.sck_count_comp.txt | cut -f1,2,3 > chr22.org_sim_count_comp.txt

# get frequency histogram of the error rate (in percetage)
cat chr22.sim_true_false_comp.txt | awk '/chr/{printf "%3.0f\n", ($4/$2)*100}' | sort -n > chr22.sim_true_false_comp.error_rate.txt
# Get frequency histogram
java -jar -Xmx1G /home/rhiea/codes/txtColumnSummary.jar 1 chr22.sim_true_false_comp.error_rate.txt | awk '{print $1" "$2}' > chr22.sim_true_false_comp.error_rate.hist


# Total error rate over origin unique kmer counts
cat chr22.sim_true_false_comp.txt | awk '/chr/{printf "%3.0f\n", ($4/$2)*100}' | sort -n > chr22.sim_true_false_comp.error_rate.txt