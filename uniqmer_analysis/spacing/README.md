# Spacing

How far apart are these uniq-mer regions?

## Dependencies
- bedtools
- canu
- R
- Java

## Inputs
- **Chromosome or Genome Name** : Must be a fasta file in the form `<Chromosome or Genome Name>.fasta`
- **Uniq-mer meryl database** : Refer to [Uniq-mer Database](https://github.com/cchu70/mashmap_postfilter/tree/master/Uniq-mer_database)
- **Chromosome gaps** : Bedfile of the gaps in the original chromosome or genome



## Command
```
/path/to/sck_spacing.sh <Chromosome> <Uniq-mer meryl DB> <Chromosome Gaps Bedfile>
```

## Output
- **Bedfile of positions of the uniq-mers** : WARNING: I hardcoded the kmer length of k=21
- **Uniq-mer Intervals** : Bedfile where overlapping kmers are merged together into a single interval
- **Gap Intervals** : Bedfiles where the uniq-mers are not found
- **Gap Intervals + Gaps in Chromosome** : Merged the gaps with existing gaps in the chromosome
- **Spacing** : Textfile listing the size of the intervals
- **Histogram Plot** : Text file to plot as a histogram. Lists the size and the frequency. 
