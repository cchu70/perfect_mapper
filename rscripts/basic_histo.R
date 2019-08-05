# Script that takes in a single column and plots a frequency histogram

library(ggplot2)
library(scales)



args = commandArgs(trailingOnly=T)

out=args[1]
table=args[2]

title=args[3]
x_header=args[4]

breaks=strtoi(args[5])

cwd = args[6]


setwd(cwd)


df = read.table(table, header="F")

# hist(df$V1, main = "True mashmap p_value frequency", xlab="p_values of true mashmap alignments", breaks=13)

hist(df$V1, main=title, xlab=x_header, breaks=breaks, ylim=c(0,ylim))