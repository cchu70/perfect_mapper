## Counting k-mers
This section takes each alignments and count the number of:
- shared uniq-mers
- non-shared uniq-mers
- shared non-uniq-mers
- non-shared non-uniq-mers

### Inputs
- **Target Fasta** : The target (or reference) that will be use to simulate nanopore reads and mapped onto by said simulated reads 
- **Simulated Reads Fasta** : simulated nanopore reads with fasta headers indicating their original start and end position in the reference
- **Alignment File** : A SAM file (no header) or mashmap output file with the **ground truth** appended at the end of each line
- **Alignment File Type** : Indicate which type. ('sam' or 'mashmap')
- **K-mer list** : A file containing each kmer that exists in the Target, or set of "true" kmers, and their corresponding frequency. There are kmers with a frequency of 0. Refer to _Unique db_ package
- **k-mer size** : Indicate the size of the k-mers used in the **k-mer list** mentioned above.

### Command
To get the k-mer counts for each alignment, run the following
```
python /path/to/weighted_jaccard/weighted_jaccard_count.py <Read Fasta> <Target Fasta> <Align File> <Align File type> <kmer list file> <k-mer size> > outfile.txt

# Example
python /path/to/weighted_jaccard/weighted_jaccard_count.py /data/Phillippy/projects/perfect-polish/chr22_info/chr22.sim_reads.fasta /data/Phillippy/projects/perfect-polish/chr22_info/chr22.fasta /data/Phillippy/projects/perfect-polish/chr22_info/representative_only.multiple_aligns_only.rev_false.minimap2_N50_30kb.real.sam /data/Phillippy/projects/perfect-polish/chr22_info/chr22.asm.sck_list.txt sam 21 > representative_only.multiple_aligns_only.rev_false.minimap2_N50_30kb.real.k_counts.0_4.txt
```

### Output
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

## Scheme Scoring
This outlines how to get the scores for the alignments given the kmer counts from the previous stage

### Inputs
- **Alignment file of Simulated Reads + Ground Truth** The ground truth ("True" or "False") is appended at the end of each line. From the **Counting K-mers** stage above. 
- **Scheme start** : Starting weight to give to the k-mers
- **Scheme end** : Last weight to give to the k-mers
- **Scheme step** : Step size between the start and end schemes. Each combination of these weights for uniq and non-uniq-mers will be tested.
- **Prefix** : for naming output files. Recommend including information on the current run.

### Commands
```
# 1.) Calculate the scores for each scheme
python /path/to/weighted_jaccard/weighted_jaccard_scheme_score.py <k-mer count alignment file> > <Scheme score file>

# Example
# python /path/to/weighted_jaccard/weighted_jaccard_scheme_score.py representative_only.multiple_aligns_only.rev_false.minimap2_N50_30kb.real.k_counts.0_4.txt 0 10 1 > representative_only.multiple_aligns_only.rev_false.minimap2_N50_30kb.real.scheme_scores.0_4.txt

# Schemes tested: (0,0), (0,1), ..., (1,0), (1,1), ...

# 2.) Evaluate the performance of each scheme
python /path/to/weighted_jaccard/weighted_jaccard_eval_main.py <Scheme score file> <prefix>

# 3.) Plot a ROC curve
Rscript /path/to/weighted_jaccard/<place holder>.R <prefix>.performance.txt
```
### Outputs
- **Scheme score file**
  1. read name
  2. Map truth ("P" or "S" for chosen mapper)
  3. Alignment start index
  4. Alignment end index
  5. Ground Truth ("True" or "False")
  6. Scores in the form `(<uniqmer_weight>,<non-uniqmer_weight>)=score`, tab delimited
```
# Example
chrX_59147642_59157542_+        S       59218397        59228062        False   0.84987067      (1.000, 22.000)=0.0633438601309 (1.000, 23.000)=0.0620884289746 (1.000, 24.000)=0.0608852041935
chrX_59147642_59157542_+        S       59783244        59793022        False   0.84966251      (1.000, 22.000)=0.0419976944886 (1.000, 23.000)=0.0410373066424
...
```
- **Incorrectly aligned reads file** `<prefix>.wrong_aligned_reads.txt`: List the reads that were incorrectly mapped
- **Performance File** `<prefix>.performance.txt`: Test file with each scheme and it's performance. Use to plot ROC curves and other performance analysis
  1. Scheme `(w_uniq,w_nonuniq)`
  2. True Positive count
  3. False Positive count
  4. False Negative count
  5. True Negative count
  6. True positive rate
  7. False positive rate
  8. Precision
 - **Plots** : Plot of the performance of each weight, labeled on each point 
 ![](images/chrX_prepolished_weighted_jaccard_whole_nums.png)
