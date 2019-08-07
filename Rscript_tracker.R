# Rscript tracking


# Script to plot the number of kmers remaining after introducing error rates
uniqmer_counts_df = read.table("uniqmer_counts.to_plot.txt")
ggplot(data=uniqmer_counts_df, aes(x = V1, y = V2, col=V3)) + geom_point() + xlab("Introduced error") + ylab("Number of uniqmers")

# To plot the retention rates
rerun_df = read.table("rerun.AAEFG_weighted_jaccard_performance.k_10.to_plot.txt")
ggplot(data=rerun_df, aes(x = V1, y = V2, col=V4)) + geom_point() + xlab("Introduced Error Rate") + ylab("Retention Rate")

# Plot retentionrates including A to B and B to A
ggplot(data=rerun_df, aes(x = V1, y = V2, col=paste(V4,V5,sep="_"))) + geom_point() + xlab("Introduced Error Rate") + ylab("Retention Rate") + theme_bw()

# Compare with minimap2

