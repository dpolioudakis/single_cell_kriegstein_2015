
# Damon Polioudakis
# 2016-04-14
# Plot bulk RNAseq of VZ from Luis and Jason's ATAC versus pooled
# scRNAseq VZ from Kriegstein 2015
################################################################################

rm(list=ls())
sessionInfo()

library(ggplot2)

### Load data and assign variables

## Load data

buExDatDF <- read.csv("../data/htseq/bulk_VZ_CP_from_ATAC/Exprs_HTSCexon.csv"
                     , row.names = 1)

scExDatDF <- read.csv("../data/htseq/Exprs_HTSCexon.csv", row.names = 1)

# Picard Sequencing Statistics - scRNAseq
picStatsScDF <- read.csv("../metadata/PicardToolsQC.csv")

# Picard Sequencing Statistics - bulk RNAseq
picStatsBuDF <- read.table("../metadata/bulk/alignment_summary.txt", fill = TRUE
                           , header = TRUE)

# GC and Gene Length
load("../../source/ENSEMBLhg19_UnionAnno.rda")
lengthGCdF <- ENSEMBLhg19.70UnionAnno
rm(ENSEMBLhg19.70UnionAnno)


## Variables

graphCodeTitle <- "Compare_Bulk_to_scRNAseq.R"
outGraph <- "../analysis/graphs/Compare_Bulk_to_scRNAseq_"
outGenesListDir <- "../analysis/GO_enrichment/Pooled_Vs_Bulk/Subsets_MAplot_0.5_-0.5_0"
dir.create(outGenesListDir, recursive = TRUE)
dir.create(file.path(outGenesListDir, "Lists"))
dir.create(file.path(outGenesListDir, "Background"))
dir.create(file.path(outGenesListDir, "Background_HighExpr"))
dir.create(file.path(outGenesListDir, "Background_LowExpr"))

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 18)))
theme_update(plot.title = element_text(size = 14))
################################################################################

### Functions


################################################################################

### Remove ERCCs
buExDatDF <- head(buExDatDF, -97)
scExDatDF <- head(scExDatDF, -97)
################################################################################

### Process data - Mean counts, Filter, read depth normalize, and pool in groups
### bulk and scRNAseq

ftScExDatDF <- scExDatDF
# ## Filter for number of mapped reads > 1.5*10^6
# # Identify Cells IDs with number of mapped reads > 1.5*10^6
# pfCells <- picStatsScDF$X[picStatsScDF$PF_READS_ALIGNED > 1.5*10^6]
# # Filter expression dataframe for cell IDs with number of mapped reads > 1.5*10^6
# ftScExDatDF <- scExDatDF[ ,colnames(scExDatDF) %in% pfCells]


## Read depth normalize

# Randomly split into 5 pooled groups

set.seed(11)
ncol(ftScExDatDF)
rNums <- sample(1:130, 130, replace = FALSE)
rndmGroups <- split(rNums, ceiling(seq_along(rNums) / (length(rNums) / 5)))
pScExDF <- data.frame(lapply(rndmGroups
                             , function (group) {apply(ftScExDatDF[ ,group], 1, sum)}))
head(pScExDF, 20)
pScIDs <- lapply(rndmGroups
                 , function (group) colnames(ftScExDatDF)[group])

# Read depth normalize pooled scRNAseq by number mapped to exons
rDep <- (apply(pScExDF, 2, sum) / 10^6)
pScTpmExonExDatDF <- pScExDF / rDep

# Read depth normalize bulk by number mapped to exons
rDep <- (apply(buExDatDF, 2, sum) / 10^6)
buTpmExonExDatDF <- buExDatDF / rDep

# Read depth normalize by number of mapped reads

# pooled scRNAseq
sel <- lapply(pScIDs, function(group) picStatsScDF$X %in% group)
splitStats <- lapply(sel, function(group) picStatsScDF[group, ])
mappedReads <- sapply(splitStats, function(group) sum(as.numeric(group$PF_READS_ALIGNED)))
pScTpmExDatDF <- pScExDF / (mappedReads / 10^6)

# Read depth normalize bulk by number mapped reads
buTpmExDatDF <- buExDatDF / (picStatsBuDF$PF_READS_ALIGNED / 10^6)


## Mean counts

# Mean counts for bulk RNAseq
mnBuEx <- apply(buTpmExDatDF, 1, mean)
lg2mnBuEx <- log(mnBuEx, 2)

# Mean counts for pooled scRNAseq groups
mnPdScEx <- apply(pScTpmExDatDF, 1, mean)
lg2mnPdScEx <- log(mnPdScEx, 2)
################################################################################

### Subset transcripts

## Use MA plot values to subset

# Added + 0.01 to all counts to prevent Inf values
maDF <- data.frame(Avg = (0.5*log((mnPdScEx + 0.01) * (mnBuEx + 0.01), 2))
                   , Log2Ratio = log(((mnPdScEx + 0.01) / (mnBuEx + 0.01)), 2))

# Split genes into groups based off expression level and bias towards bulk
# or scRNAseq
biasDF <- rbind(
  data.frame(
  # Low expressed bias towards bulk
  # Transcripts > 0.5 M (log ratio) towards bulk and > 0 A (mean average)
  # (expression > ~1 TPM)
    Genes = maDF[maDF$Avg < 0 & maDF$Log2Ratio < -0.5, ]
    , Gene_Set = "Low Expression, Higher Bulk")
  
  , data.frame(
  # Low expressed bias towards scRNAseq
  # Transcripts > 0.5 M (log ratio) towards pooled scRNAseq and > 0 A (mean
  # average) (expression > ~1 TPM)
    Genes = maDF[maDF$Avg < 0 & maDF$Log2Ratio > 0.5, ]
    , Gene_Set = "Low Expression, Higher Pooled")
  
  , data.frame(
  # High expressed bias towards bulk
  # Transcripts > 0.5 M (log ratio) towards bulk and > 0 A (mean average)
  # (expression > ~1 TPM)
    Genes = maDF[maDF$Avg > 0 & maDF$Log2Ratio < -0.5, ]
    , Gene_Set = "Higher Expression, Higher Bulk")
  
  , data.frame(
  # How expressed bias towards scRNAseq
  # Transcripts > 0.5 M (log ratio) towards pooled scRNAseq and > 0 A (mean
  # average) (expression > ~1 TPM)
    Genes = maDF[maDF$Avg > 0 & maDF$Log2Ratio > 0.5, ]
    , Gene_Set = "Higher Expression, Higher Pooled")
)
################################################################################

### Plot Length and GC content for gene subsets

# Add length and gc content
biasDF <- merge(biasDF, lengthGCdF, by.x = "row.names", by.y = "row.names")

# Calc means - to add to plot
meansLength <- aggregate(UnionExonLength ~ Gene_Set, biasDF, mean)
meansGC <- aggregate(UnionGCcontent ~ Gene_Set, biasDF, mean)

# Calc number of samples - to add to plot
nSamples <- data.frame(table(biasDF$Gene_Set))
colnames(nSamples) <- c("Gene_Set", "Number_of_Samples")


## Barplot number of samples in each subset

ggplot(nSamples, aes(y = Number_of_Samples, x = Gene_Set, fill = Gene_Set)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none") +
  xlab("Gene Sets") +
  ylab("Number of Samples") +
  ggtitle(paste0("Compare_Bulk_to_scRNAseq.R"
               , "\nNumber of Genes in Subsets with Biased Expression"
               , "\n"))
ggsave(paste0(outGraph, "Subsets_Numbers.pdf"))


## Length Bias Plot

ggplot(biasDF, aes(y = UnionExonLength, x = Gene_Set, col = Gene_Set)) +
  geom_boxplot(outlier.shape = NA) +
  # Adjust limits after outlier removal
  coord_cartesian(ylim = range(boxplot(biasDF[
    biasDF$Gene_Set == "Higher Expression, Higher Bulk", ]$UnionExonLength
    , plot = FALSE)$stats) * c(.9, 1.1)) +
  geom_text(data = meansLength, aes(label = paste("Mean:", round(UnionExonLength, 2)), y = UnionExonLength)) +
  geom_text(data = nSamples, aes(label = paste("n =", Number_of_Samples)
                                 , y = meansLength$UnionExonLength - 200)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none") +
  xlab("Gene Sets") +
  ylab("Length") +
  ggtitle(paste0("Compare_Bulk_to_scRNAseq.R"
    , "\nGene Length for Subsets of Genes with Biased Expression"
    , "\nP-value: Low Expression, Higher Bulk vs Higher Pooled scRNA-seq: "
      ,signif(t.test(x = biasDF[biasDF$Gene_Set == "Low Expression, Higher Bulk", ]$UnionExonLength
      , y = biasDF[biasDF$Gene_Set == "Low Expression, Higher Pooled", ]$UnionExonLength)$p.value, 3)
    ,"\nP-value: High Expression, Higher Bulk vs Higher Pooled scRNA-seq: "
      ,signif(t.test(x = biasDF[biasDF$Gene_Set == "Higher Expression, Higher Bulk", ]$UnionExonLength
      , y = biasDF[biasDF$Gene_Set == "Higher Expression, Higher Pooled", ]$UnionExonLength)$p.value, 3)
    , "\n")
)
ggsave(paste0(outGraph, "Subsets_Length.pdf"))


## GC Bias Plot

# Filter genes with NA GC content (2 genes)
ggDF <- biasDF[! is.na(biasDF$UnionGCcontent), ]

ggplot(ggDF, aes(y = UnionGCcontent, x = Gene_Set, col = Gene_Set)) +
  geom_boxplot(outlier.shape = NA) +
  # Adjust limits after outlier removal
  coord_cartesian(ylim = range(boxplot(ggDF[
    ggDF$Gene_Set == "Higher Expression, Higher Bulk", ]$UnionGCcontent
    , plot = FALSE)$stats) * c(.9, 1.1)) +
  geom_text(data = meansGC, aes(label = paste("Mean:", round(UnionGCcontent, 2)), y = UnionGCcontent + 0.02)) +
  geom_text(data = nSamples, aes(label = paste("n =", Number_of_Samples)
                                 , y = meansGC$UnionGCcontent - 0.02)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none") +
  xlab("Gene Sets") +
  ylab("GC Content") +
  ggtitle(paste0("Compare_Bulk_to_scRNAseq.R"
                 , "\nGC Content for Subsets of Genes with Biased Expression"
                 , "\nP-value: Low Expression, Higher Bulk vs Higher Pooled scRNA-seq: "
                 ,signif(t.test(x = ggDF[ggDF$Gene_Set == "Low Expression, Higher Bulk", ]$UnionGCcontent
                                , y = ggDF[ggDF$Gene_Set == "Low Expression, Higher Pooled", ]$UnionGCcontent)$p.value, 3)
                 ,"\nP-value: High Expression, Higher Bulk vs Higher Pooled scRNA-seq: "
                 ,signif(t.test(x = ggDF[ggDF$Gene_Set == "Higher Expression, Higher Bulk", ]$UnionGCcontent
                                , y = ggDF[ggDF$Gene_Set == "Higher Expression, Higher Pooled", ]$UnionGCcontent)$p.value, 3)
                 , "\n")
)
ggsave(paste0(outGraph, "Subsets_GCcontent.pdf"))
################################################################################

### GO Analysis of subsets of genes

# Format gene names to Ensembl
fixEnsNamesDF <- biasDF
fixEnsNamesDF$Row.names <- gsub("\\..*", "", fixEnsNamesDF$Row.names)

# Split dataframe by gene subset
dfList <- split(fixEnsNamesDF, fixEnsNamesDF$Gene_Set)


## Background gene lists

# All genes
# Format background gene list names
bgListAll <- gsub("\\..*", "", names(mnPdScEx))

# High expressed
# > 0 A (mean average) (expression > ~1 TPM)
bgListHigh <- row.names(maDF)[maDF$Avg > 0]
# Format background gene list names
bgListHigh <- gsub("\\..*", "", names(mnPdScEx))

# Low expressed
# < 0 A (mean average) (expression > ~1 TPM)
bgListLow <- row.names(maDF)[maDF$Avg < 0]
# Format background gene list names
bgListLow <- gsub("\\..*", "", names(mnPdScEx))


## Write out tables of genes for each subset

write.table(data.frame(
  ensembl_gene_id = dfList$"Low Expression, Higher Bulk"$Row.names
    , systemCode = "En")
  , file = paste(outGenesListDir, "/lists/LowExpr_HigherBulk.txt", sep = "")
  , row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
)
write.table(data.frame(
  ensembl_gene_id = dfList$"Low Expression, Higher Pooled"$Row.names, systemCode = "En")
            , file = paste(outGenesListDir, "/lists/LowExpr_HigherPooled.txt", sep = "")
            , row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
)
write.table(data.frame(
  ensembl_gene_id = dfList$"Higher Expression, Higher Bulk"$Row.names
    , systemCode = "En")
  , file = paste(outGenesListDir, "/lists/HighExpr_HigherBulk.txt", sep = "")
  , row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
)
write.table(data.frame(
  ensembl_gene_id = dfList$"Higher Expression, Higher Pooled"$Row.names
    , systemCode = "En")
  , file = paste(outGenesListDir, "/lists/HighExpr_HigherPooled.txt", sep = "")
  , row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
)

# Write out table of all genes (for background)
write.table(data.frame(ensembl_gene_id = bgListAll, systemCode = "En")
            , file = paste(outGenesListDir, "/Background/background.txt", sep = "")
            , row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
)
# Write out table of high expressed genes (for background)
write.table(data.frame(ensembl_gene_id = bgListHigh, systemCode = "En")
            , file = paste(outGenesListDir, "/Background_HighExpr/background.txt", sep = "")
            , row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
)
# Write out table of low expressed genes (for background)
write.table(data.frame(ensembl_gene_id = bgListLow, systemCode = "En")
            , file = paste(outGenesListDir, "/Background_LowExpr/background.txt", sep = "")
            , row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
)


## RUN GO ELITE ##

## Plotting the GO Output

Make_GO_Elite_File_List <- function (filesDir) {
  goEliteFiles <- list.files(filesDir)
  idx <- grep(".*GO_z-score_elite.txt", goEliteFiles)
  goEliteFiles <- goEliteFiles[idx]
  head(read.csv(file = paste0(goElitePath, goEliteFiles[1]), sep = "\t"))$Ontology.Type
  goEliteFiles
}

Plot_GO_Elite_Results <- function (goEliteFileNamesL) {
  for(i in 1:length(goEliteFileNamesL)){
    thismod = goEliteFileNamesL[i]
    nameThisMod = goEliteFileNamesL[i]
    print(nameThisMod)
    tmp = read.csv(file = paste0(goElitePath, goEliteFileNamesL[i]), sep = "\t")
    # Loop through ontology types
    lapply(unique(tmp$Ontology.Type), function(GOtype) {
      tmp = subset(tmp, Ontology.Type==GOtype)
      tmp = tmp[ ,c(2,4,9)] ## Select GO-terms and Z-score
      tmp = tmp[order(tmp$Z.Score, decreasing = TRUE), ] #
      tmp = tmp[which(tmp[ ,2] >= 1), ]
      
      if (nrow(tmp) == 0){
        print(paste("No enriched GO-terms for:", nameThisMod))
        next
      } else if (nrow(tmp) < 10) {
        tmp1 = tmp ## Take top 10 Z-score
        tmp1 = tmp1[order(tmp1$Z.Score),] ##Re-arrange by increasing Z-score
        par(mar=c(4,22,4,4))
        barplot(tmp1$Z.Score, horiz = TRUE, col = "blue"
                , names.arg = tmp1$Ontology.Name, cex.names=1.2, las=1
                , main = paste("Gene Ontology Plot of", GOtype
                               , "\n", thismod)
                , xlab = "Z-Score")
        abline(v=2,col="red")
        next
      } else {
        tmp1 = tmp[c(1:20), ] ## Take top 10 Z-score
      }
      tmp1 = tmp1[order(tmp1$Z.Score), ] ##Re-arrange by increasing Z-score
      par(mar = c(4, 22, 4, 4))
      barplot(tmp1$Z.Score, horiz = TRUE, col = "blue"
              , names.arg = tmp1$Ontology.Name
              , cex.names = 1.2, las = 1
              , main = paste("Gene Ontology Plot of", GOtype
                             , "\n", thismod)
              , xlab = "Z-Score")
      abline(v = 2, col = "red")
      cat('Done ...', thismod, GOtype, '\n')
    })
  }
}

# Using all genes as background
goEliteFileList <- Make_GO_Elite_File_List("../analysis/GO_enrichment/Subsets_MAplot_0.5_-0.5_0/GO_Elite_Output/GO-Elite_results/CompleteResults/ORA_pruned/")
pdf("../analysis/graphs/GO_Elite.pdf", height = 6, width = 8)
Plot_GO_Elite_Results(goEliteFileList)
dev.off()

# Using high expressed genes as background for high expresse
goEliteFileList <- Make_GO_Elite_File_List("../analysis/GO_enrichment/Subsets_MAplot_0.5_-0.5_0/GO_Elite_Output_HighExpr/GO-Elite_results/CompleteResults/ORA_pruned/")
pdf("../analysis/graphs/GO_Elite_HighExpr.pdf", height = 6, width = 8)
Plot_GO_Elite_Results(goEliteFileList)
dev.off()

# Using low expressed genes as background for low expressed
goEliteFileList <- Make_GO_Elite_File_List("../analysis/GO_enrichment/Subsets_MAplot_0.5_-0.5_0/GO_Elite_Output_LowExpr/GO-Elite_results/CompleteResults/ORA_pruned/")
pdf("../analysis/graphs/GO_Elite_LowExpr.pdf", height = 6, width = 8)
Plot_GO_Elite_Results(goEliteFileList)
dev.off()






uniquemodcolors = unique(colors)
uniquemodcolors = uniquemodcolors[-which(uniquemodcolors=="grey")]

pdf(outGOgraphs, height = 6, width = 8)

for(i in 1:length(uniquemodcolors)){
  thismod= uniquemodcolors[i]
  nameThisMod = paste("MM",uniquemodcolors[i],sep="")
  print(nameThisMod)
  tmp=read.csv(file = paste(outModsLists
                            , "/GO-Elite_results/CompleteResults/ORA/"
                            , nameThisMod, "-GO_z-score_elite.txt", sep = "")
               , sep = "\t")
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