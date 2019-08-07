# Rscript tracking


# Script to plot the number of kmers remaining after introducing error rates
uniqmer_counts_df = read.table("uniqmer_counts.to_plot.txt")
ggplot(data=uniqmer_counts_df, aes(x = V1, y = V2, col=V3)) + geom_point() + xlab("Introduced error") + ylab("Number of uniqmers")

