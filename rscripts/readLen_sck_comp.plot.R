# Plot dot plot of readLen vs unique kmer counts, colorized based on
# which sequence came from which group

library(ggplot2)
library(scales)

setwd("/data/chucl/chr22_tests/")

args = commandArgs(trailingOnly=T)

out=args[1]
table=args[2]
error=strtoi(args[3])
x_size=strtoi(args[4])
y_size=strtoi(args[5])


readLen_sck=read.table(table, header=F)

ggplot(data=readLen_sck, aes(x=V2, y=V3, color=V1)) + geom_point() + geom_abline(intercept = 0, slope=error)

ggsave(file=paste(out,'.scatter_plot.png', sep=""), height = x_size, width = y_size)