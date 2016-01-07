# Damon Polioudakis
# 2015-12-04
# Graphs of Kriegstein 2015 single-cell RNA-seq alignment statistics
rm(list=ls())
sessionInfo()

library(ggplot2)
library(reshape2)

# Input gene expression and metadata
exDatDF <- read.csv("../data/mmc3.csv", header = TRUE, row.names = 1)
metDatDF <- read.csv("../data/mmc2.csv", header = TRUE)

# Number of genes expressed
gExprd <- apply(exDatDF, 2, function(genes) length(subset(genes, genes > 0)))
gExprd <- sort(gExprd)
mgExprd <- mean(gExprd)
pdf("../analysis/graphs/Alignment_Stats_Graph_Genes_Expressed.pdf")
barplot(gExprd, xlab = "Cells", ylab = "Number of genes TPM > 0", xaxt = 'n'
        , main = paste("Kriegstein 2015: Number of Genes Expressed"
                       ,"\nMean: ", mgExprd, sep=""))
Axis(side = 1, labels = FALSE)
dev.off()

# Number of cells in which gene expression > CPM - stat Kriegstein used in paper
nCellExprD <- apply(exDatDF, 1, function(genes) length(subset(genes, genes > 1)))
nCellExprD <- sort(nCellExprD)
nkriegFilt <- sum(nCellExprD > 2)
tGene <- nrow(exDatDF)
mnCellExprD <- mean(nCellExprD)
pdf("../analysis/graphs/Alignment_Stats_Graph_Genes_Expr_Per_Cell.pdf")
barplot(nCellExprD, xlab = "Genes", ylab = "Number of cells in which gene > 1 CPM", xaxt = 'n'
        , main = paste("Kriegstein 2015: Number of Cells With Gene Expressed"
                       ,"\nMean number of cells in which gene > 1 CPM: ", mnCellExprD
                       ,"\nTotal Genes: ", tGene
                       ,"\nKriegstein Filter (Red Line: Genes > 1 CPM in > 2 Cell):", nkriegFilt, sep=""))
abline(a = 2, b = 0, col = "red")
Axis(side = 1, labels = FALSE)
dev.off()

# Library size and read pairs aligned
pdf("../analysis/graphs/Alignment_Stats_Graph_Library_Size_Alignment.pdf")
nReadsMapDF <- melt(metDatDF[c("Cell", "LibrarySize..Read.Pairs..Million.", "MappedPairs")])
mMapped <- mean(metDatDF$MappedPairs)
mLibSize <- mean(metDatDF$LibrarySize..Read.Pairs..Million.)
ggplot(nReadsMapDF, aes(x = Cell, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_discrete(breaks = NULL) +
  xlab("Cells") +
  ylab("Read Pairs (Million)") +
  ggtitle(paste("Kriegstein 2015: Library Size and Mapped Pairs"
                ,"\nMean Library Size: ", mLibSize
                ,"\nMean Mapped Pairs: ", mMapped, sep = ""))
dev.off()

# Percent aligned
pdf("../analysis/graphs/Alignment_Stats_Graph_Percent_Aligned.pdf")
alignRate <- sort(as.numeric(substr(metDatDF$AlignmentRate..Pairs,1,4)))
mAlignRate <- mean(alignRate)
barplot(alignRate, xlab = "Cells", ylab = "Percent Aligned (Pairs)"
        , main = paste("Kriegstein 2015: Alignment Rate"
                       ,"\nMean Alignment Rate: ", mAlignRate, sep=""))
dev.off()
