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
rerun_minimap2_df = read.table("AAEFG.outfiles.performance.to_plot.txt")
> ggplot(data=rerun_minimap2_df, aes(x = V1, y = V2, col=paste(V4,V5,sep="_"))) + geom_point() + xlab("Introduced Error Rate") + ylab("Retention Rate") + theme_bw()
> 
> rerun_minimap2_df$type = "minimap2"
> rerun_df$type = "WJ"
> 
> rerun_combined <- rbind(rerun_minimap2_df, rerun_df)
> ggplot(data=rerun_combined, aes(x = V1, y = V2, col=paste(V4,V5,type,sep="_"))) + geom_point() + xlab("Introduced Error Rate") + ylab("Retention Rate") + theme_bw()
> ggplot(data=rerun_combined, aes(x = V1, y = V2, col=type)) + geom_point() + xlab("Introduced Error Rate") + ylab("Retention Rate") + theme_bw()
> 


require(gridExtra)
Loading required package: gridExtra
> grid.arrange(minimap2_perf, wj_0)
> grid.arrange(minimap2_perf, wj_0, ncol=2)
> wj_0 <- ggplot(data=subset(vary_weights, (V6 == 0)), aes(x=V1, y=V2, col=paste(V3,V4))) + geom_point() + ylim(0,1) + theme_bw() + labs(fill = "Weight 0 for non-uniq-mers")
> minimap2_perf <- ggplot(data=rerun_minimap2_df, aes(x=V1, y=V2, col=paste(V4,V5))) + geom_point() + ylim(0,1) + theme_bw() + labs(fill="Minimap2")
> grid.arrange(minimap2_perf, wj_0, ncol=2)
> wj_1 <- ggplot(data=subset(vary_weights, (V6 == 1)), aes(x=V1, y=V2, col=paste(V3,V4))) + geom_point() + ylim(0,1) + theme_bw() + labs(fill = "Weight 1 for non-uniq-mers")
> grid.arrange(minimap2_perf, wj_0, wj_1, ncol=3)
> wj_0 <- ggplot(data=subset(vary_weights, (V6 == 0, V5 == 2)), aes(x=V1, y=V2, col=paste(V3,V4))) + geom_point() + ylim(0,1) + theme_bw() + labs(fill = "Weight 0 for non-uniq-mers")
Error: unexpected ',' in "wj_0 <- ggplot(data=subset(vary_weights, (V6 == 0,"
> wj_0 <- ggplot(data=subset(vary_weights, (V6 == 0 & V5 == 2)), aes(x=V1, y=V2, col=paste(V3,V4))) + geom_point() + ylim(0,1) + theme_bw() + labs(fill = "Weight 0 for non-uniq-mers")
> ggplot(data=subset(vary_weights, (V6 == 0 & V5 == 2)), aes(x=V1, y=V2, col=paste(V3,V4))) + geom_point() + ylim(0,1) + theme_bw() + labs(fill = "Weight 0 for non-uniq-mers")
> wj_0_2 <- ggplot(data=subset(vary_weights, (V6 == 0 & V5 == 2)), aes(x=V1, y=V2, col=paste(V3,V4))) + geom_point() + ylim(0,1) + theme_bw() + labs(fill = "Weight 0 for non-uniq-mers")
> wj_0_4 <- ggplot(data=subset(vary_weights, (V6 == 0 & V5 == 4)), aes(x=V1, y=V2, col=paste(V3,V4))) + geom_point() + ylim(0,1) + theme_bw() + labs(fill = "Weight 0 for non-uniq-mers")
> wj_0_8 <- ggplot(data=subset(vary_weights, (V6 == 0 & V5 == 8)), aes(x=V1, y=V2, col=paste(V3,V4))) + geom_point() + ylim(0,1) + theme_bw() + labs(fill = "Weight 0 for non-uniq-mers")
> wj_0_16 <- ggplot(data=subset(vary_weights, (V6 == 0 & V5 == 16)), aes(x=V1, y=V2, col=paste(V3,V4))) + geom_point() + ylim(0,1) + theme_bw() + labs(fill = "Weight 0 for non-uniq-mers")
> wj_0_32 <- ggplot(data=subset(vary_weights, (V6 == 0 & V5 == 32)), aes(x=V1, y=V2, col=paste(V3,V4))) + geom_point() + ylim(0,1) + theme_bw() + labs(fill = "Weight 0 for non-uniq-mers")
> grid.arrange(wj_0_2, wj_0_4, wj_0_8, wj_0_16, wj_0_32, ncol = 3)
> ggplot(data=subset(vary_weights, (V6 == 0 & V5 == 16)), aes(x=V1, y=V2, col=paste(V3,V4))) + geom_point() + ylim(0,1) + theme_bw() + labs(fill = "Weight 0 for non-uniq-mers")
> 
> 
> weight_1 = read.table("weights_1.AAEFG_weighted_jaccard_performance.txt")
> ggplot(data = weight_1, aes(x=V1, y=V2, col=paste(V5,V6))+ geom_point() + ylim(0,1) + theme_bw() 
+ 

> ggplot(data = weight_1, aes(x=V1, y=V2, col=paste(V5,V6)) + geom_point() + ylim(0,1) + theme_bw() 
+ 

> ggplot(data = weight_1, aes(x=V1, y=V2, col=paste(V5,V6))) + geom_point() + ylim(0,1) + theme_bw() 
Warning message:
Removed 2640 rows containing missing values (geom_point). 
> weight_1 = read.table("weights_1.AAEFG_weighted_jaccard_performance.txt")
> ggplot(data = weight_1, aes(x=V1, y=V2, col=paste(V5,V6))) + geom_point() + ylim(0,1) + theme_bw() 


