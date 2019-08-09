# Read length versus Uniq-mer Count

How does length of simulated nanopore reads impact the number of uniq-mers remaining

## Dependencies

## Steps
1. Make a uniq-mer + position meryl database (refer to [Uniq-mer database](https://github.com/cchu70/mashmap_postfilter/tree/master/Uniq-mer_database). 
2. Dump the reads on to the uniq-mer database
```
meryl-lookup -dump -sequence <reads> -mers <uniq-mer database> | awk '$3 == "T" {print $0}' > <dump>
```
3. For each read, count how many of it's k-mers are unique. 
```
cat <dump> | awk 'BEGIN{ curr_read = "", count = 0 } { if($1 != curr_read) { print $curr_read"\t"count; curr_read = $1, count = 0 } else { count = count + 1} } > <Read + uniq-mer count>
```
4. From here I'd imagine doing some more fancy awk statements to parse the read names to get the start and end positions of the read to get the length but I don't really have time to test it at the moment
