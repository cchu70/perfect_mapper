# R script to plot do plots with linear regression lines


library(ggplot2)
library(scales)



args = commandArgs(trailingOnly=T)

out=args[1]
table=args[2]

title=args[3]
x_header=args[4]
y_header=args[5]

cwd = args[6]

header=args[7]


setwd(cwd)

isHeaderTrue <- function(x) {
	if (x == "T") {return(T)}
	else {return(F)}
}

ggplotRegression <- function (fit, xl, yl) {
	require(ggplot2)
	ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
		ggtitle(title) +
		xlab(xl) + ylab(yl) +
		geom_point() +
		stat_smooth(method = "lm", col = "red") +
		labs(title = paste(paste(title, "\n", sep =""), "Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
			"Intercept =",signif(fit$coef[[1]],5 ),
			" Slope =",signif(fit$coef[[2]], 5),
			" P =",signif(summary(fit)$coef[2,4], 5)))
}


df = read.table(table, header=isHeaderTrue(header))

fit <- lm(df[,2] ~ df[,1], data = df)

ggplotRegression(fit, x_header, y_header)

ggsave(file=paste(out,'.basic_dotplot.png', sep=""))