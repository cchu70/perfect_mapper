# Rscript to plot the performance of weighted jaccard weighting schemes, table output from Simulated Reads Weighted Jaccard Performance Percentages

library(ggplot2)

args = commandArgs(trailingOnly=T)


out=args[1]				# Output name (no type ending)
table=args[2]			# Weighted Jaccard Performance Percentage output


minimap2_perf=read.table(table, header=T)

# Combine the combination of aligned reads and the part errored 
minimap2_perf$Combination = paste(minimap2_perf$V4, minimap2_perf$V5)

# Plot
ggplot(data=minimap2_perf, aes(x=V1, y=V2, col=Combination)) + geom_point() + ylim(0,1) + xlab("Error rate introduced") + ylab("Retention Rate") + theme_bw()
# Subset or change the coloring schemes to view different combinations

ggsave(file=paste(out,'.plot_sim_error_minimap2.png', sep=""))