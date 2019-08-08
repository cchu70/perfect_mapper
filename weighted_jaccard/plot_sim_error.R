# Rscript to plot the performance of weighted jaccard weighting schemes, table output from Simulated Reads Weighted Jaccard Performance Percentages

library(ggplot2)

args = commandArgs(trailingOnly=T)


out=args[1]				# Output name (no type ending)
table=args[2]			# Weighted Jaccard Performance Percentage output


false_uniq_rates=read.table(table, header=T)

ggplot(data=false_uniq_rates, aes(x=V2) + geom_histogram() + xlab(x_label) + ylab(y_label)

ggsave(file=paste(out,'.freq_hist.png', sep=""))