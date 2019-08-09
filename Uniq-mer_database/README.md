# Collecting the uniq-mers
We collect uniq-mers by using canu's meryl package, which allows quick identification of uniq-mers.

# Dependencies
- canu (meryl, meryl-lookup)

# Make the uniq-mer databases
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
