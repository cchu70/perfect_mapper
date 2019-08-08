# Weighted Jaccard Testing
This collection of scripts calculates and evaluates the performance of a weighted Jaccard Scheme to select the correct alignment.

![Weighted Uniq-mer Jaccard Index](images/NIH_SIP_Poster_Images-Weighting.png)

# Dependencies
- canu
- bedtools
- minimap2 or mashmap

# Using an existing reference
This section consists of the following steps:
1. Parsing an alignment file (BAM or mashmap output)
2. Evaluating mapper
3. Evaluating performance of multiple weighting schemes
4. Comparing performance

# Ground Truth
To determine if an alignment is true, we treat the alignment of the non-errored version of a simulated read as the Ground Truth, and as long as the simulated read aligns such that it covers at least 50% of the original alignment, this counts as true

## Inputs
- **Target Fasta** : Fasta file used to derive the reads
- **Unique meryl database** : Database with only unique kmers and their positions in the reference mentioned above (refer to make meryl db script)
- **Simulated Read Alignment File** : A Sam file with the alignments for the simulated reads onto the target sequence
- **Origin Read Fasta**
- **Prefix** : To name files

## Commands

1. Align the origin reads onto the target sequence
```
/path/to/weighted_jaccard/get_new_ground_truth.sh <origin reference fasta> <origin read fasta> <prefix> 
```
2. Label the alignments of the simulated reads alignments with the ground truth file
```
samtools view <simulated read alignment file> | /path/to/ground_truth_from_origin_alignment.py <Ground Truth File> > output
```
## Outputs
- **Alignment file for origin reads on the target**
- **Ground Truth File**
  1. Read name
  2. Index of start of the alignment
- **Alignment file of Simulated Reads + Ground Truth**

# Something Else

## Inputs
- **Target Fasta** : The target (or reference) that will be use to simulate nanopore reads and mapped onto by said simulated reads 
- **Simulated Reads Fasta** : simulated nanopore reads with fasta headers indicating their original start and end position in the reference
- **Simulated Reads True Positions bedfile** : Bedfile with the original positions of each simulated read
```
# Origin_seq    start   end   read_name
chr22   2000    8000    chr22_2000_8000_+
```
## Outputs

# Counting k-mers
This section takes each alignments and count the number of:
- shared uniq-mers
- non-shared uniq-mers
- shared non-uniq-mers
- non-shared non-uniq-mers

## Inputs
- **Target Fasta** : The target (or reference) that will be use to simulate nanopore reads and mapped onto by said simulated reads 
- **Simulated Reads Fasta** : simulated nanopore reads with fasta headers indicating their original start and end position in the reference
- **Alignment File** : A SAM file (no header) or mashmap output file with the **ground truth** appended at the end of each line
- **Alignment File Type** : Indicate which type. ('sam' or 'mashmap')
- **K-mer list** : A file containing each kmer that exists in the Target, or set of "true" kmers, and their corresponding frequency. There are kmers with a frequency of 0. Refer to _Unique db_ package
- **k-mer size** : Indicate the size of the k-mers used in the **k-mer list** mentioned above.

## Command
To get the k-mer counts for each alignment, run the following
```
python /path/to/weighted_jaccard/weighted_jaccard_count.py <Read Fasta> <Target Fasta> <Align File> <Align File type> <kmer list file> <k-mer size> > outfile.txt

# Example
python /path/to/weighted_jaccard/weighted_jaccard_count.py /data/Phillippy/projects/perfect-polish/chr22_info/chr22.sim_reads.fasta /data/Phillippy/projects/perfect-polish/chr22_info/chr22.fasta /data/Phillippy/projects/perfect-polish/chr22_info/representative_only.multiple_aligns_only.rev_false.minimap2_N50_30kb.real.sam /data/Phillippy/projects/perfect-polish/chr22_info/chr22.asm.sck_list.txt sam 21 > representative_only.multiple_aligns_only.rev_false.minimap2_N50_30kb.real.k_counts.0_4.txt
```

## Output
- **Count File** : A text file with reduced information from the original alignment file and the k-mer counts
 1. Read Name
 2. "P" - Primary, "S" - Secondary
 3. Start alignment index in the target sequence
 4. End alignment index in the target sequence
 5. Ground Truth
 6. Shared uniq-mers count
 7. Shared non-uniq-mers count
 8. Non-shared uniq-mers count
 9. Non-shared non-uniq-mers count
 10. Shared erroneous k-mers count (Not shown in example, was a later addition)
 11. Non-shared erroneous k-mers count (Not shown in example, was a later addition)
 
```
# representative_only.multiple_aligns_only.rev_false.minimap2_N50_30kb.real.k_counts.0_4.txt
chr22_part11_139296_176509      P       12166854        12204055        True    802     9445    2223    50752
chr22_part11_139296_176509      S       11788202        11825524        False   6       919     3947    67519
chr22_part02_8424_15157 P       10843069        10849804        True    24      2132    91      9016
chr22_part02_8424_15157 S       15872490        15879225        False   0       2118    134     9049
chr22_part24_103255_113647      P       18586770        18597158        True    42      2971    183     14026
chr22_part24_103255_113647      S       21526317        21536678        False   0       2863    195     14322
...
```

# Simulation on the GAGE locus
To test the effects of variable error rates, we simulated random error in one half of the GAGE locus, and then realigned simulated reads back onto this new GAGE locus.

![](images/NIH_SIP_Poster_Images-Simulated_error_test.png)


## Inputs

## Outputs
