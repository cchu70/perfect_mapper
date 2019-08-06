# Rscript to plot the output file from false_uniq_analysis.sh ${prefix}.uniqmer_loss.to_plot.txt

library(ggplot2)

args = commandArgs(trailingOnly=T)


out=args[1]				# Output name (no type ending)
table=args[2]			# ${prefix}.uniqmer_loss.to_plot.txt
x_label=args[3]			# X-axis label
y_label=args[4]			# Y-axis label


uniq_comp=read.table(table, header=T)

ggplot(data=uniq_comp, aes(x=V2, y=V3)) + geom_point(size=0.2) + geom_abline(intercept = 0, slope=1) + xlab(x_label) + ylab(y_label)

ggsave(file=paste(out,'.dot_plot.png', sep=""))