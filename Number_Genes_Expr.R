# Damon Polioudakis
# 2016-01-03
# Graph of number of genes expressed in 2015-12 scRNAseq run C196-001
# SxaQSEQsXbp060L2 and 2015-12 scRNAseq additional sequencing C196-002
rm(list=ls())
sessionInfo()

# Input gene expression table (counts from HTseq)
exDatDF <- read.csv("../data/HTSC/Exprs_HTSCexon.csv", header = TRUE, row.names = 1)
dim(exDatDF)

# Number of genes expressed
gExprd <- apply(exDatDF, 2, function(genes) length(subset(genes, genes >= 1)))
mgExprd <- mean(gExprd)
mean(gExprd)
median(gExprd)
pdf("../analysis/graphs/Histogram_Genes_Expressed.pdf")
png("../analysis/graphs/Histogram_Genes_Expressed.png")
hist(gExprd, breaks = 20, xlab = "Genes Detected with >=1 reads"
	, main = paste("Number of Genes Expressed", "\nMean: ", mgExprd, sep=""))
dev.off()


