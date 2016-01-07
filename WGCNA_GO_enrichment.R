sessionInfo()
rm(list=ls())

library(WGCNA)
library(biomaRt)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

################################################################################

# Compare modules to meta data

load("../data/WGCNA_3_Modules_SP12.rda")

metDatDF <- read.csv("../data/mmc2.csv", header = TRUE)
# Each column of modulesColors is modules made with different parameters
modulesColors <- modulesColors[ ,7]


################################################################################

# merged = mergeCloseModules(exprData= exDatDF, colors = geneTree$labels, cutHeight=dthresh)
colors = vector(mode="list")
# colors = cbind(colors, labels2colors(modulesColors))
colors = modulesColors
colors = as.character(colors)
## Find module membership for genes
MEs0 = moduleEigengenes(exDatDF, modulesColors)$eigengenes
MEs = orderMEs(MEs0)
modNames = substr(names(MEs), 3, length(names(MEs)))
geneModuleMembership = as.data.frame(cor(exDatDF, MEs, use = "p"))
nSamples = dim(exDatDF)[1]
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")



# Assemble module .txt's,background .txt, and make lists of most correlated
# genes with ME in modules

# First assemble hubs list, for each module a data frame with most significantly
# correlated genes at the top with p values
hubs=list()
nModules <- length(unique(modulesColors))
for(i in 1:nModules){
  
  hubs[[i]]=list(data = as.data.frame(sort(geneModuleMembership[,i],decreasing=T,index.return=TRUE)))
  names = rownames(geneModuleMembership)[hubs[[i]]$data[,2]]
  hubs[[i]]$data = as.data.frame(hubs[[i]]$data[,-2])
  rownames(hubs[[i]]$data) = names
  
  idx = match(rownames(hubs[[i]]$data),rownames(MMPvalue))
  hubs[[i]]$data = cbind(hubs[[i]]$data,MMPvalue[idx,i])  
  colnames(hubs[[i]]$data) = c(names(geneModuleMembership)[i],names(MMPvalue)[i])
}

names(hubs) = names(geneModuleMembership)

# Now background list

background = colnames(exDatDF)

# Now module lists

modules = list()
for(i in 1:nModules){
  modules[[i]] = subset(background,colors==modNames[i])
  names(modules)[i] = paste("MM",modNames[i],sep="")
}

## Now GOElite

AddEnsembl <- function (geneList) {
  moduleGenes <- data.frame(geneList)
  # bioMart manual:
  # http://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
  # Attribute: values to retrieve
  # Filters: input query type
  # Values: input query
  #ensembl <- useMart("ensembl")
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
  ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
  # Data frame of module Ensembl IDs and gene symbols
  moduleEnsemblDF <- getBM(  attributes = c("ensembl_gene_id", "hgnc_symbol")
                             , filters = "hgnc_symbol"
                             , values = moduleGenes
                             , mart = ensembl
  )
  moduleEnsemblDF
}
background <- background <- AddEnsembl(background)
modules <- lapply(modules[1:nModules-1], function(x) AddEnsembl(x))


for (i in 1:(nModules-1)){
  write.table(data.frame(sourceID = modules[[i]][1], systemCode = "En")
              , file=paste("../analysis/GO_enrichment/lists/", names(modules)[i],".txt", sep="")
              , row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
}

write.table(data.frame(sourceID = background[1], systemCode = "En"),file="../analysis/GO_enrichment/background/background.txt",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")




####Plotting the GO Output

uniquemodcolors = unique(colors)
uniquemodcolors=uniquemodcolors[-which(uniquemodcolors=="grey")]

pdf("../analysis/graphs/GOElite_BP_plot_Modules.pdf",height=8,width=8)

for(i in 1:length(uniquemodcolors)){
  thismod= uniquemodcolors[i]
  nameThisMod = paste("MM",uniquemodcolors[i],sep="")
  print(nameThisMod)
  tmp=read.csv(file=paste("../analysis/GO_enrichment/output/GO-Elite_results/CompleteResults/ORA_pruned/",nameThisMod,"-GO_z-score_elite.txt",sep=""),sep="\t")
  tmp=subset(tmp,Ontology.Type=="biological_process")
  tmp=tmp[,c(2,4,9)] ## Select GO-terms and Z-score
  tmp=tmp[order(tmp$Z.Score,decreasing=T),] #
  tmp=tmp[which(tmp[,2] >= 1),]
  
  if (nrow(tmp) == 0){
    print(paste("No enriched GO-terms for:", nameThisMod))
    next
  } else if (nrow(tmp)<10) {
    tmp1=tmp ## Take top 10 Z-score
    tmp1 = tmp1[order(tmp1$Z.Score),] ##Re-arrange by increasing Z-score
    par(mar=c(4,22,4,4))
    barplot(tmp1$Z.Score,horiz=T,col="blue",names.arg= tmp1$Ontology.Name,cex.names=1.2,las=1,main=paste("Gene Ontology Plot of",thismod,"Module"),xlab="Z-Score")
    abline(v=2,col="red")
    next
  } else {
    tmp1=tmp[c(1:10),] ## Take top 10 Z-score
  }
  tmp1 = tmp1[order(tmp1$Z.Score),] ##Re-arrange by increasing Z-score
  par(mar=c(4,22,4,4))
  barplot(tmp1$Z.Score,horiz=T,col="blue",names.arg= tmp1$Ontology.Name,cex.names=1.2,las=1,main=paste("Gene Ontology Plot of",thismod,"Module"),xlab="Z-Score")
  abline(v=2,col="red")
  
  cat('Done ...',thismod,'\n')
}

dev.off()



tmp=read.csv(file=paste("../analysis/GO_enrichment/output/GO-Elite_results/CompleteResults/ORA_pruned/MMblue-GO_z-score_elite.txt",sep=""),sep="\t")
tmp=subset(tmp,Ontology.Type=="biological_process")
tmp=tmp[,c(2,4,9)] ## Select GO-terms and Z-score
tmp=tmp[order(tmp$Z.Score,decreasing=T),] #
tmp=tmp[which(tmp[,2] >= 10),]
