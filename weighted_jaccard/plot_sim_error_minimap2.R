# Rscript to plot the performance of weighted jaccard weighting schemes, table output from Simulated Reads Weighted Jaccard Performance Percentages

library(ggplot2)

args = commandArgs(trailingOnly=T)


out=args[1]				# Output name (no type ending)
table=args[2]			# Weighted Jaccard Performance Percentage output


minimap2_perf=read.table(table, header=T)

minimap2_perf$Combination = paste(minimap2_perf$V4, minimap2_perf$V5)

# Plot
ggplot(data=minimap2_perf, aes(x=V1, y=V2, col=Combination)) + geom_point() + ylim(0,1) + xlab("Error rate introduced") + ylab("Retention Rate") + theme_bw()
# Subset or change the coloring schemes to view different combinations

# Examples
# Plot only the performances of weighting schemes where the non-uniqmer weight is 0
# ggplot(data=subset(vary_weights, (V6 == 0)), aes(x=V1, y=V2, col=paste(V3,V4))) + geom_point() + ylim(0,1) + theme_bw() + labs(fill = "Weight 0 for non-uniq-mers")

ggsave(file=paste(out,'.plot_sim_error_minimap2.png', sep=""))