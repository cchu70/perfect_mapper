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
