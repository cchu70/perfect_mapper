# Collecting the uniq-mers
We collect uniq-mers by using canu's meryl package, which allows quick identification of uniq-mers.

# Dependencies
- canu (meryl, meryl-lookup)

# Quick run: Make the uniq-mer databases
## Inputs
-  **Prefix** : Usually the chromosome number, but just used to label the output files
- **Input Fasta** : Fasta file to count k-mers. Usually used to get the k-mers of a polished genome or full chromosome
- **K-mer size** : Length of each k-mer. I used k=21
- **Threads** : Number of threads to run meryl
- **Memory** : Amount of memory to run meryl

## Command
```
/path/to/make_uniq_meryl_db.sh <prefix> <input fasta> <kmer size> <threads> <memory>

# Example
./make_uniq_meryl_db.sh chrX chrX.polished.fasta 21 40 60g
```

## Outputs
- **Main meryl database** : Meryl database containing all the kmers in the input fasta and their corresponding values (counts)
- **Equal-to 1 meryl database** : Meryl database only containing the kmers in the main meryl DB that has a value (count) of 1
- **Uniq-mer positions** : Text file containing each k-mer in the Equal-to-1 meryl DB and it's position index in the input fasta
- **Uniq-mer position meryl database** : Meryl database containing the uniq-mers and their corresponding value (position index)
- **K-mer list** : Text file containing each k-mer in the input fasta file and it's corresponding count, tab delimited. This file is an input to some scripts in the `weighted_jaccard` folder. 


# Meryl Databases Tutorial
The following lists how to make meryl databases

## Uniq-mer database of some sequence

1. Count the kmers
```
meryl count k=$k_size <input_fa> output XXX.main.meryl
```
2. Filter only the kmers that occur once (uniq-mers)
```
meryl equal-to 1 XXX.main.meryl output XXX.equal_to_1.meryl
```

## Uniq-mer + Position database
1. Take the output of **Uniq-mer database of some sequence** and dump the original fasta file used to make this database to get the position of each kmer
```
meryl-lookup -dump -mers XXX.equal_to_1.meryl -sequence <input_fa> -threads <threads> -memory <memory> | awk -F "\t" '$5==1{print $4"\t"$2}'  > XXX.uniq_pos.to_import.txt
```
2. Import into a new database
```
meryl-import -k <k-mer size> -kmers XXX.uniq_pos.to_import.txt -output XXX.uniq_pos.meryl -threads $threads -memory $memory
```


## Uniq-mer database for a set of reads
This is useful to reduce the set of uniq-mers to consider 
1. Count kmers in the reads
```
meryl count k=21 <XXX.reads.fasta> output XXX.reads.meryl
```
2. Make unique position meryl database of the target sequence. Refer above to **Uniq-mer database + position**
3. Intersect the read database with the target
```
meryl intersect XXX.uniq_pos.meryl XXX.reads.meryl output XXX.reads.intersect.meryl
```



