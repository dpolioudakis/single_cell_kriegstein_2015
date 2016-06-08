# Damon Polioudakis
# 2016-01-03
# Graph of number of genes expressed in 2015-12 scRNAseq run C196-001
# SxaQSEQsXbp060L2 and 2015-12 scRNAseq additional sequencing C196-002
rm(list=ls())
sessionInfo()

# Input gene expression table (counts from HTseq)
exDatDF <- read.csv("../data/htseq/Exprs_HTSCexon.csv", header = TRUE, row.names = 1)
dim(exDatDF)

# Picard Sequencing Statistics - scRNAseq
picStatsDF <- read.csv("../metadata/PicardToolsQC.csv")

# Calculate CPM
exDatDF <- exDatDF[ ,order(names(exDatDF))]
cpmDF <- exDatDF / (picStatsDF$PF_READS_ALIGNED / 10^6)

# Number of genes expressed - raw counts
gExprd <- apply(exDatDF, 2, function(genes) length(subset(genes, genes >= 1)))
mgExprd <- mean(gExprd)
mean(gExprd)
median(gExprd)
hist(gExprd, breaks = 20, xlab = "Genes Detected with >=1 reads"
	, main = paste("Number of Genes Expressed", "\nMean: ", mgExprd, sep=""))
dev.off()

# Number of genes expressed - CPM
gExprd <- apply(cpmDF, 2, function(genes) length(subset(genes, genes >= 1)))
mgExprd <- mean(gExprd)
mean(gExprd)
median(gExprd)
pdf("../analysis/graphs/Number_Genes_Detected_Histogram.pdf")
hist(gExprd, breaks = 20, col = "blue", xlab = "Number of genes with >=1 CPM"
	, main = paste("Number of Genes Detected Per Cell", "\nMean: ", round(mgExprd,2), sep=""))
dev.off()