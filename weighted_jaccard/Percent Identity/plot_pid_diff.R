# Rscript to plot the percent Identity Differences between true alignments and primary alignments



library(ggplot2)

args = commandArgs(trailingOnly=T)


out=args[1]				# Output name (no type ending)
table=args[2]			# pid_correctness file (minimap2 correctness or weighted jaccard correctness file)


pid_diff=read.table(table, header=T)

ggplot(data = pid_diff, aes(x = V1, y=V2)) + geom_jitter(position=position_jitter(w=0.4), alpha=0.5, color="darkblue", size=1.5) + xlab("Correct Primary Alignment") + ylab("%idy difference") + theme_bw()

ggsave(file=paste(out,'.plot_pid_diff.png', sep=""))