# Damon Polioudakis
# 2016-06-07
# Number of reads per cell
################################################################################

rm(list=ls())

library(xlsx)

# Load data and assign variables

# RNA STAR Stats Lane 1
metDatDF <- read.xlsx2("../metadata/Cell paper - updated attributes.xlsx", 1)
################################################################################

# Combine total reads stats from Lane 1 and Lane 2
totReads <- metDatDF$LibrarySize..Read.Pairs..Million.
# Convert from factor to numeric
totReads <- as.numeric(levels(totReads))[totReads]
totReads <- totReads * 10^6

mnReads <- mean(totReads)

# Plot as histogram
pdf("../analysis/graphs/Number_Reads_Per_Cell_Histogram.pdf")
hist(totReads, breaks = 20, col = "blue", xlab = "Number of reads per cell"
     , main = paste("Number of Reads Per Cell", "\nMean: ", round(mnReads,1), sep=""))
dev.off()