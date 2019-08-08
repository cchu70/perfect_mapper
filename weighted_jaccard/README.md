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

## Ground Truth
To determine if an alignment is true, we treat the alignment of the non-errored version of a simulated read as the Ground Truth, and as long as the simulated read aligns such that it covers at least 50% of the original alignment, this counts as true

### Inputs
- **Target Fasta** : Fasta file used to derive the reads
- **Unique meryl database** : Database with only unique kmers and their positions in the reference mentioned above (refer to make meryl db script)
- **Simulated Read Alignment File** : A Sam file with the alignments for the simulated reads onto the target sequence
- **Origin Read Fasta**
- **Prefix** : To name files

### Commands
```
# 1. Align the origin reads onto the target sequence

/path/to/weighted_jaccard/get_new_ground_truth.sh <origin reference fasta> <origin read fasta> <prefix> 

# 2. Label the alignments of the simulated reads alignments with the ground truth file

samtools view <simulated read alignment file> | /path/to/ground_truth_from_origin_alignment.py <Ground Truth File> > output
```
### Outputs
- **Alignment file for origin reads on the target**
- **Ground Truth File**
  1. Read name
  2. Index of start of the alignment
- **Alignment file of Simulated Reads + Ground Truth** The ground truth ("True" or "False") is appended at the end of each line

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

## Percent Identity Effects
To see the potential effects of variable percent identity between a read and multiple candidate alignments, I calculate the percent ID for each alignment. Then for each mapping method, I compared the percent ID between the mapper's "primary" alignment and the true alignment (if said primary alignment was false), or the second best alignment (if the primary alignment was already true to begin with) 

Optional, but you can just analyze reads with multiple alignments. Refer to the **Format** folder for more info. 

### Inputs
- **SAM/BAM file + ground truth** : Alignment file to be evaluated. Will use the cigar string to calculate the pid. To append the ground truth on each alignment line, refer to the **Ground Truth** section. 
- **Reduced SAM/BAM file + ground truth + scheme score** : Run the **Scheme Score** section. Select only one weighting scheme to test. Select the best scheme by looking at the performance curves, mentioned in the same section

### Commands
```
# 1.) Append the %idy on each alignment string
samtools view [filter options] <SAM or BAM file> | python /path/to/weighted_jaccard/calc_percent_identity.py > <Reduced SAM + pid>


# 2.) Evaluate the performance of minimap2
python /path/to/weighted_jaccard/pid_diff_eval/py <SAM + pid> > <minimap2 correctness>


# 3.) Evaluate the performance of the weighted jaccard method using selected scheme

python /path/to/weighted_jaccard/weighted_jaccard_pid_diff_eval_schemes.py <Reduced SAM/BAM file + ground truth + scheme score> > <Weighted Jaccard correctness>

# 4.) Plot
Rscript /path/to/plot_pid_diff.R <Weighted Jaccard correctness Prefix for png file> <Weighted Jaccard correctness>
Rscript /path/to/plot_pid_diff.R <minimap2 correctness prefix for png file> <minimap2 correctness>

```
### Outputs
- **Reduced SAM + pid** : Reduced alignment file with the percent identity appended to the end
  1. Read Name
  2. Mapper truth ("P" for primary, "S" for secondary)
  3. Alignment start index
  4. Alignment end index
  5. Ground Truth
  6. Percent Identity
- **minimap2 correctness** and **Weighted Jaccard correctness** : For each read the mapper aligned, lists if it is correct, and it's corresponding percent identity difference between either the True alignment (if it exists) or the second best alignment if it is already true
- **minimap2 correctness plot** and **Weighted Jaccard correctness plot** : Plots of the difference in percent identity of multiple candidate alignments

![](images/chrX_minimap2_pid_diff.png)
![](images/chrX_wj_pid_diff_plot.png)


## Simulation on the GAGE locus
To test the effects of variable error rates, we simulated random error in one half of the GAGE locus, and then realigned simulated reads back onto this new GAGE locus.

![](images/NIH_SIP_Poster_Images-Simulated_error_test.png)

### Inputs
Steps 1 - 2 : Evaluate Minimap2
- **GAGE_A reference fasta** : Fasta file to part A of GAGE or of target sequence
- **GAGE_B reference fasta** : Fasta file to part B of GAGE or of target sequence
- **GAGE_A Simulated reads fasta** : Simulated reads originating from GAGE_A reference sequence
- **GAGE_B Simulated reads fasta** : Simulated reads originating from GAGE_B reference sequence
- **Error Rate Start** : Starting error rate to introduce
- **Error Rate End** : Ending error rate to introduce
- **Error Rate Step** : Increments of error rates to test between the Start and End error rates
- **Iterations** : Number of times to test at each error rate
- **Prefix** : The run name to track files
- **Script path** : Path to the script `test_alignment_on_sim_errors.sh`

Step 3 : Evaluate Weighted Jaccard
- **Errors String** : String of the error rates to evaluate, comma delimited (no spaces). Ensure that the errors listed exist as directories in the form `error_<error rate>/`. 
- **Prefix** : The run you want to evaluate (ie pick the one you used above)
- **Uniqmer Weights list** : String comma delimited (no spaces) list of the weights to try on uniqmers
- **Non-uniqmer Weights list** : String comma delimited (no spaces) list of the weights to try on non-uniqmers

### Commands
```
# 1.) Simulate errors and create alignments

python /path/to/weighted_jaccard/simulate_error_main.py <GAGE_A> <GAGE_B> <GAGE_A_reads> <GAGE_B_reads> <Error start> <Error end> <Error step> <Iterations> <prefix> /path/to/test_alignment_on_sim_errors.sh

  # Example
  # python /path/to/weighted_jaccard/simulate_error_main.py GAGE_A.fasta GAGE_B.fasta GAGE_A.sim_reads.fasta GAGE_B.sim_reads.fasta 0.0 0.15 0.1 10 XXX /path/to/test_alignment_on_sim_errors.sh

# 1_alt.) To run the simulation on a single condition (like only introduce errors to GAGE_A at only 0.1 error rate)

/path/to/test_alignment_on_sim_errors.sh <GAGE_A> <GAGE_B> <Which to error ("A" or "B")> <Which NOT to error ("A" or "B")> <Error rate> <Iterations> <GAGE_A_reads> <GAGE_B_reads> <prefix> 

  # Example
  /path/to/test_alignment_on_sim_errors.sh GAGE_A.fasta GAGE_B.fasta A B 0.1 5 GAGE_A.sim_reads.fasta GAGE_B.sim_reads.fasta XXX

# 2.) Evaluate the output files from minimap2

python /path/to/sim_error_eval.py <File of file paths to all the out files> > <minimap2 Performance Percentage>.to_plot

# 3.) Evaluate the weighted Jaccard approach on the alignment files produced from Step 1

python /path/to/weighted_jaccard_sam_input_k_count.py <Errors String> <Prefix> <Uniq weights> <non-uniq weights> > <Weighted Jaccard Performance Count>

  # Example
  python /path/to/weighted_jaccard_sam_input_k_count.py 0,0.1,0.2,0.3,0.4,0.5,0.6 XXX 1,2,4,8 0,1 > weighted_jaccard_performance.txt
  
# 4.) Plot

cat <Weighted Jaccard Performance Count> | grep "GAGE" | awk '{print $1"\t"($2/($2 + $3))"\t"$5"\t"$6"\t"$7"\t"$8}' > <Weighted Jaccard Performance Percentage>.to_plot.txt


Rscript /path/to/plot_sim_error_weighted_jaccard.R <Weighted Jaccard Performance Percentage>.to_plot.txt
Rscript /path/to/plot_sim_error_minimap2.R <minimap2 Performance Percentage>.to_plot.txt

```

### Outputs
Step 1:
- **Error rate directories** : Each error rate will get it's own directory in the form `error_<error rate>/`. For each iteration `${counter}`:
- **Errored fasta file** : Fasta file of the part that we introduce error to. In the form `${prefix}.err_${error_rate}_${which_to_error}.v_${counter}.fasta`
- **Split Target fasta file** : Concatenated target fasta with the fasta record of the Part with no error and the fasta record of the Errored Part fasta file. In the form `${prefix}_split.err_${error_rate}_${which_to_error}.v_${counter}.fasta`
- **Part A alignment file** : minimap2 alignment of reads from A onto the **Split Target fasta file** with parameters `-N50 -r3000`. In the form `${prefix}_minimap2.N50_r3k.split.err_${error_rate}_${which_to_error}.v_${counter}.aligned_A.sam`
- **Part B alignment file** : minimap2 alignment of reads from B onto the **Split Target fasta file** with parameters `-N50 -r3000`. In the form `${prefix}_minimap2.N50_r3k.split.err_${error_rate}_${which_to_error}.v_${counter}.aligned_B.sam`
- **Performance Out file** : Documents the how many reads from what part mapped where for the corresponding error rate. In the form `${prefix}.err_${error_rate}_${which_to_error}.out`
  1. Error rate
  2. Which Part Errored 'A/B'
  3. Number of reads from A aligned to Part A
  4. Number of reads from A aligned to Part B
  5. Number of reads from B aligned to Part B
  6. Number of reads from B aligned to Part A
  7. Number of read from A unaligned
  8. Number of reads from B unaligned
  9. Part A alignment file name : For tracking purposes
  10. Part B alignment file name : For tracking purposes 

Step 2:
- **minimap2 performance percentage** : Text file documenting the performance for each error rate
  1. Error rate
  2. Retention Rate: Percentage of reads from column 5 mapping to the correct part (ex. Percentage of reads from Part A mapping back to Part A)
  3. Which part was Errored ("GAGE_A")
  4. Reads from which part was aligned
  
Step 3:
- **Weighted Jaccard Performance** : Text file documenting the performance of each weighting scheme for each error rate for each iteration
  1. Error rate
  2. Correct count : The number of reads from column 6 that were mapped to the correct part
  3. Incorrect count : The number of reads from column 6 that were mapped to the incorrect part
  4. Unmapped count : The number of reads without an alignment
  5. Which part was errored : GAGE_A or GAGE_B
  6. Reads from which part was aligned : GAGE_A or GAGE_B
  7. Unique Kmer weight
  8. Non-unique kmer weight
  9. Corresponding Sam file (For tracing)

Step 4:
- **Weighted Jaccard Performance Percentage** : Text file containing information on the performance on each error rate for each weighting scheme
  1. Error rate
  2. Retention Rate (0 - 1): Number of reads in column 4 mapped to the correct part
  3. Which part errored
  4. Reads from which part aligned
  5. Uniq-mer weight
  6. Non-uniq-mer weight
- **Plots!** : In the form `*.plot_sim_error_minimap2.png` or `*.plot_sim_error_weighted_jaccard.png`. Enough information are in the above plots to compare the two, view specific weights, etc. 

![](images/GAGE_vary_weights_performance.plot_sim_error_weighted_jaccard.png)


