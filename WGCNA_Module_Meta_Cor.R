sessionInfo()
rm(list=ls())

library(WGCNA)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

################################################################################

# Compare modules to meta data

load("../data/WGCNA_3_Modules.rda")

metDatDF <- read.csv("../data/mmc2.csv", header = TRUE)

# Define numbers of genes and samples
nGenes = ncol(exDatDF);
nSamples = nrow(exDatDF);
# Recalculate MEs with color labels
# Each column of modulesColors is modules made with different parameters
MEs0 = moduleEigengenes(exDatDF, modulesColors[ ,1])$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, metDatDF, use = "p");

metDatDF$AlignmentRate..Pairs <- as.numeric(
  substr(as.character(metDatDF$AlignmentRate..Pairs), 1, 4))
rSqM <- matrix(NA, nrow = ncol(MEs), ncol = (ncol(metDatDF) - 1))
for(ME in 1:ncol(MEs)) {
  for(metDat in 2:ncol(metDatDF)) {
    rSqM[ME, metDat-1] <- summary(lm(MEs[ ,ME] ~ metDatDF[ ,metDat]))$adj.r.squared
  }
}
moduleTraitCor <- sign(rSqM) * abs(rSqM)^(1/2)

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(metDatDF)[-1],
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))













# Select specific traits

gene.sigs=matrix(NA,nrow=7,ncol=ncol(exDatDF)) # create a vector to hold the data
for(i in 1:ncol(gene.sigs)) {
  
  exprvec= as.numeric(t(exDatDF)[,i]) # get the expression vector for ith gene
  
  anatSrc = sqrt(max(summary(lm(exprvec~as.factor(metDatDF[ ,2])))$adj.r.squared,0))
  age = bicor(as.numeric(metDatDF[ ,3]), exprvec)# calculate r correlation value for numeric variables
  totReads = bicor(metDatDF[ ,4], exprvec)# calculate r correlation value for numeric variables
  numMap = bicor(metDatDF[ ,5], exprvec)# calculate r correlation value for numeric variables
  alignRate = bicor(as.numeric(metDatDF[ ,6]), exprvec)# calculate r correlation value for numeric variables
  genesTagged = bicor(metDatDF[ ,6], exprvec)# calculate r correlation value for numeric variables
  # Well IDs are numbers, ie: "160535191 160535175 160091869 160091634"
  
  gene.sigs[, i]=c(anatSrc, age, totReads, numMap, alignRate, genesTagged)
}


gene.sigs[1,] =numbers2colors(as.numeric(gene.sigs[1,]),signed=FALSE,centered=FALSE,blueWhiteRed(100),lim=c(0,1)) # For categorical
gene.sigs[2,] =numbers2colors(as.numeric(gene.sigs[2,]),signed=FALSE,centered=FALSE,blueWhiteRed(100),lim=c(0,1)) # For categorical
gene.sigs[3,] =numbers2colors(as.numeric(gene.sigs[3,]),signed=FALSE,centered=FALSE,blueWhiteRed(100),lim=c(0,1)) # For categorical
gene.sigs[4,] =numbers2colors(as.numeric(gene.sigs[4,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
rownames(gene.sigs)=c("anatSrc","age","totReads", "numMap")

# Parameters used to construct module [,3]: DS= 2, MMS= 30, MEcor= 0.25
modules.traits=data.frame(modules.colors[,3]
                          , gene.sigs[1,]
                          , gene.sigs[2,]
                          , gene.sigs[3,]
                          , gene.sigs[4,])
modules.traits.labels=c(module.labels[3],rownames(gene.sigs))

pdf("../analysis/dendro_modules_traits_corr.pdf",height=25,width=20)
plotDendroAndColors(gene.tree,modules.traits,groupLabels=modules.traits.labels,addGuide=TRUE,dendroLabels=FALSE,main="Dendro and traits correlation")
dev.off()
