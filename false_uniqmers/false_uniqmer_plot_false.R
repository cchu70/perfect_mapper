# Rscript to plot the frequency of rate of false uniqmers, table output from false_uniq_analysis.sh in the form ${prefix}.false_uniqmer_rate.to_plot.txt

library(ggplot2)

args = commandArgs(trailingOnly=T)


out=args[1]				# Output name (no type ending)
table=args[2]			# ${prefix}.uniqmer_loss.to_plot.txt
x_label=args[3]			# X-axis label
y_label=args[4]			# Y-axis label


false_uniq_rates=read.table(table, header=T)

ggplot(data=false_uniq_rates, aes(x=V2) + geom_histogram() + xlab(x_label) + ylab(y_label)

ggsave(file=paste(out,'.freq_hist.png', sep=""))