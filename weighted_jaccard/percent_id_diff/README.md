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
python /path/to/weighted_jaccard/pid_diff_eval.py <SAM + pid> > <minimap2 correctness>


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

![](../images/chrX_minimap2_pid_diff.png)
![](../images/chrX_wj_pid_diff_plot.png)
