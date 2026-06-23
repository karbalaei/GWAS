# 5. Results visualization 
Results visualization are crucial to explain and inspect the results. You can look at your results using Manhattan plot and QQ-plot. Also it is good to check the Lambda value from your GWAS.


## Manhattan plot

See here a great way of making a Manhattan plot:
https://github.com/ipdgc/Manhattan-Plotter

## QQ-plot

Making the QQ plot is a way to vizualize the potential inflation

```
data = read.table("YOUR_summaray_stats.TBL",header=T)
# optional subsetting e.g. number of datasets or MAF
# newdata <- subset(data, HetDf > 9)
# newdata <- subset(data, MAF > 0.05)
observed <- sort(p)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
# can update the name of output file if needed here, can also change to pdf("file.pdf")
png("qqplot.png")
# note that the range is here set to 10 on both X and Y axis you might want to change this if you have more significant hits
plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
dev.off()
```

## Lambda value
```
data = read.table("YOUR_summaray_stats.TBL",header=T)
# optional subsetting e.g. number of datasets or MAF
# newdata <- subset(data, HetDf > 9)
# newdata <- subset(data, MAF > 0.05)
p <- newdata$Pvalue
n <- length(newdata$Pvalue)
x2obs <- qchisq(as.numeric(as.character(p)), 1, lower.tail = FALSE)
x2exp <- qchisq((1:n)/n, 1, lower.tail = FALSE)
lambda <- median(x2obs)/median(x2exp)
lambda
# lambda1000 can be used to normalize unbalanced case-control numbers
# change the Ncases and Ncontrols
lambda1000 <- 1 + ( lambda -1 ) * (1/'Ncases' + 1/'Ncontrols')/(1/1000 + 1/1000)
lambda1000
```











