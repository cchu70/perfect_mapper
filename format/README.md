# Formatting
Ways to trim down alignment files to filter out certain alignments depending on what you want to look at.

For example, I only want to consider the reads' alignments if for each of those reads, there are at least two candidate alignments to pick from


## Merge bam alignments
WARNING: Probably deprecated. Output of some parsing functions may not be compatible anymore.
This is if you want to represent split alignments under single alignment line. 

### Inputs
- **SAM/BAM file** : Alignment file, unsorted. This is so supplementary reads immediately precede their representative alignment

### Commands
```
samtools view [options] <SAM/BAM> | python /path/to/merge_bam_alignments.py > <Reduced Alignment file>
```
### Outputs
- **Reduced Alignment file** : Containing merged split alignments 
  1. Read name
  2. Alignment type ("P" or "S")
  3. Start index of alignment
  4. End index of alignment
  5. Mapping Quality of the representative read

## Get Multiple alignments only
Just filter out the alignments for reads that have at least 2 candidate alignments. 

### Inputs
- **Alignment file** : Any alignment file where the read name is listed in the first column

### Command
```
python /path/to/get_multiple_alignments_only.py <Alignment file> > <Multiple alignments only>
```

### Output
- **Multiple alignments only** : Same format as the input file, but only lists the alignments for reads that have at least two candidate alignments


## Get Multiple alignments only + Ground Truth

This both prints out only alignments for reads that have multiple alignments, with the ground truth appended at the end

### Inputs
- **SAM of mashmap out file** : Alignment file of mapping simulated reads to the target sequence to filter. Assumes that the alignment file is in order such that all the alignments of a given read are adjacent
- **True position bed file** : Bedfile of the true start and end positions of the simulated reads in the target reference. Assumes only one record in the target sequence

### Commands
```
# If using SAM file

samtools view [options] <SAM file> | python /path/to/get_multiple_alignments_from_bam.py <True Position bedfile> > <Multiple Alignments + Ground Truth>

# If using Mashmap output file
python /path/to/get_multiple_alignments_from_mashmap.py <mashmap out> <True Position bedfile> > <Multiple Alignments + Ground Truth>
```
### Outputs
- **Multiple Alignments + Ground Truth**:
  1. The original alignment line
  2. Ground Truth ("True" or "False")

