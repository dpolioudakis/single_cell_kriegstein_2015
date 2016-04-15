# Damon Polioudakis
# 2016-03-09
# Calculate gene length and GC content for exon union of each gene in list
# If GC content and exon union have already been calculated for reference,
# can skip to last section and load
################################################################################

rm(list = ls())
sessionInfo()

source("http://bioconductor.org/biocLite.R")
# biocLite("Genominator")
library(Repitools)
library(BSgenome)
# Run if BSgenome.Hsapiens.UCSC.hg19 is not installed:
# biocLite("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
# Run if Genominator is not installed:
# biocLite("Genominator")
library(Genominator)
library(biomaRt)
library(ggplot2)

options(stringsAsFactors = FALSE)

## Load data and assign variables

exDatDF <- read.csv("../data/htseq/Exprs_HTSCexon.csv")

# Get Gencode 18 gtf file - this was cleaned by selecting only the columns
# containing the word "exon" and with the relevant information - chrnum,
# feature type, strand, start, end, and ensg and ense IDs separated by a semicolon
gtfInfoDF <- read.table("../../source/gencode.v19.annotation.gtf", sep = "\t")
################################################################################

### Format and Filter

# Keep only the exon level features
keep <- gtfInfoDF[ , 3]=="exon"
gtfInfoDF <- gtfInfoDF[keep, ]

# Split the semicolon separated information
geneExonInfo <- unlist(strsplit(gtfInfoDF[ , 9], "[;]"))

# Finding which has gene_id for ensembl ID
genCol<- which(regexpr("gene_id ", geneExonInfo) > 0)
getSeq <- geneExonInfo[genCol]
ENSGID <- substr(getSeq, 9, 100)
length(unique(ENSGID)) #57912

transCol <- which(regexpr("transcript_id ", geneExonInfo) > 0)
tranSeq <- geneExonInfo[transCol]
ENSEID <- substr(tranSeq, 16, 100)
length(unique(ENSEID)) #196612
gtfInfoDF <- cbind(gtfInfoDF[ , c(1:8)], ENSGID, ENSEID)

geneCol <- which(regexpr("gene_name ", geneExonInfo) > 0)
geneSeq <- geneExonInfo[geneCol]
GENEID <- substr(geneSeq, 12, 100)
length(unique(GENEID)) #55763

gtfInfoDF <- cbind(gtfInfoDF[ , c(1:8)], ENSGID, ENSEID)
# 6 and 8 columns are blank - remove
gtfInfoDF <- gtfInfoDF[ , -c(6, 8)]

######## From Viveks script - I think this incorrectly only keeps 1 exon per gene
# ## Keep only one copy of each ENSEID - the gtf file records one copy for each transcript id
# keep <- match(unique(ENSEID),ENSEID)
# gtfInfo1dF <- gtfInfoDF[keep,]
# ##gtfInfoDF[,1] <- substr(gtfInfoDF[,1],4,10) ## 672406 exons is exactly what biomaRt contains
gtfInfo1dF <- gtfInfoDF
################################################################################

### Recode things for the Genominator package

# Using as.factor to coerce chromosome names can really botch things up.. beware!
# So go ahead and convert MT, X, and Y to numbers throughout, unless necessary
# for other purposes
chrnums <- gtfInfo1dF[ ,1]
chrnums[chrnums=="MT"] <- "25"
chrnums[chrnums=="X"] <- "23"
chrnums[chrnums=="Y"] <- "24"
# If there are no chr to remove this does not work
## removing Non-annotated(NT) chromosomes
# rmChR.col1 <- which(regexpr("HG", chrnums) > 0)
# rmChR.col2 <- which(regexpr("GL", chrnums) > 0) 
# rmChR.col3 <- which(regexpr("HS", chrnums) > 0)
# rmChR.col <- c(rmChR.col1, rmChR.col2, rmChR.col3)
# gtfInfo1dF <- gtfInfo1dF[-rmChR.col, ]
# chrnums <- chrnums[-rmChR.col]
gtfInfo1dF[ ,1] <- chrnums ## Check here

gtfInfoDF <- gtfInfo1dF

strinfo <- gtfInfoDF[ ,6]
strinfo[strinfo=="+"] <- 1L
strinfo[strinfo=="-"] <- -1L
gtfInfoDF[ ,6] <- strinfo

# chr integer, strand integer (-1L,0L,1L), start integer, end integer, ensg
# and transcript id
geneDat1 <- gtfInfoDF[ ,c(1 ,6 ,4 ,5 ,7 ,8)]
geneDat1 <- data.frame(as.numeric(chrnums)
                       , as.numeric(geneDat1[ ,2])
                       , as.numeric(geneDat1[ ,3])
                       , as.numeric(geneDat1[ ,4])
                       , geneDat1[ ,5]
                       , geneDat1[ ,6])
names(geneDat1) <- c("chr","strand","start","end","ensembl_gene_id","ensembl_exon_id")

geneDatX <- geneDat1[order(geneDat1[ ,1], geneDat1[ ,3]), ]
# Remove NAs from ERCC chromosomes
geneDatX <- geneDatX[complete.cases(geneDatX), ]
# Have genominator check if this is a valid data object
validAnnotation(geneDatX)
# Should take a few minutes !!!!
geneDatX <- makeGeneRepresentation(annoData = geneDat1
                                   , type = "Ugene"
                                   , gene.id = "ensembl_gene_id"
                                   , transcript.id = "ensembl_exon_id"
                                   , verbose = TRUE)

save(geneDatX, file = "../../source/Genominator_Union_Exon_Models_ENSEMBLhg19.rda")
load(file = "../../source/Genominator_Union_Exon_Models_ENSEMBLhg19.rda")
################################################################################

### Now use the genominator output to calculate GC content - Use mac laptop

geneDat2 <- cbind(geneDatX, geneDatX[ , 3] - geneDatX[ , 2])
geneDat2 <- geneDat2[order(geneDat2[ , 5]) , ]

# Change formatting again
chrNums <- geneDat2[  ,"chr"]
chrNums[chrNums=="25"] <- "M" ## important as UCSC codes MT as M
chrNums[chrNums=="23"] <- "X"
chrNums[chrNums=="24"] <- "Y"
stInfo <- geneDat2[ ,"strand"]
stInfo[stInfo==c(-1)] <- "-"
stInfo[stInfo==c(1)] <- "+"

# Calculate GC content using the union exon ranges determined by
# selecting "exons" from the gtf file above
# Convert to a genomic ranges object
gcQuery <- GRanges(paste("chr", chrNums, sep = "")
                   , IRanges(geneDat2[ , 2], geneDat2[ , 3]), strand = stInfo)
gcContent <- gcContentCalc(x = gcQuery, organism = Hsapiens)

# Take a length weighted average of GC content percentages to get the GC content
# for the union gene model
head(geneDat2)
geneDat2 <- cbind(geneDat2, gcContent)
geneDat2 <- cbind(geneDat2, gcContent * geneDat2[ , 6])
unionGenes <- by(geneDat2[ , 6], as.factor(geneDat2[ , 5]), sum)
unionGC <- by(geneDat2[ , 8], as.factor(geneDat2[ , 5]), sum)
geneDat3 <- cbind(unionGenes, unionGC / unionGenes)
colnames(geneDat3) <- c("UnionExonLength","UnionGCcontent")
ENSEMBLhg19.70UnionAnno <- geneDat3

## Save for further usage
save(ENSEMBLhg19.70UnionAnno, file="../../source/ENSEMBLhg19_Exon_Union_Anno.rda")
load("../../source/ENSEMBLhg19_Exon_Union_Anno.rda")
################################################################################

### Check Lengths and GC content to longest isoform for each gene id from biomart

# Pull Length and % GC content from Biomart, includes UTRs and CDS
ensemblMart <- useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")
ensemblMart <- useDataset("hsapiens_gene_ensembl", mart = ensemblMart)
martDF <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id"
                           , "transcript_length", "percentage_gc_content")
                , mart = ensemblMart)

# Select longest isoform for each gene id
sel <- ave(martDF$transcript_length, martDF$ensembl_gene_id
                 , FUN = max) == martDF$transcript_length
martMaxDF <- martDF[sel, ]

# Merge biomart dataframe with union exon dataframe
df <- ENSEMBLhg19.70UnionAnno
row.names(df) <- gsub("\\..*", "", row.names(df))
df <- merge(df, martMaxDF, by.x = "row.names", by.y = "ensembl_gene_id")

## Plot biomart longest isoform length and gc content vs union exon

# Length
ggplot(df, aes(x = UnionExonLength, y = df$transcript_length)) +
  geom_point(alpha = 0.25, shape = 1) +
  xlab("Union Exon Gene Length") +
  ylab("Biomart Longest Isoform Gene Length") +
  ggtitle(paste0("Calc_Union_Exon_Length_and_GC.R"
                 ,"\nCompare Union Exon Gene Length to Biomart Longest Isoform"
                 ,"\nPearson:", round(cor(df$UnionExonLength
                                          , df$transcript_length
                                    , method = "pearson"), 2)
                 ,"\nSpearman:", round(cor(df$UnionExonLength
                                           , df$transcript_length
                             , method = "spearman"), 2)))
ggsave("../analysis/Calc_Union_Exon_Length_and_GC_Compare_Length_To_Biomart.pdf")

# GC
ggplot(df, aes(x = UnionGCcontent, y = df$percentage_gc_content/100)) +
  geom_point(alpha = 0.25, shape = 1) +
  xlab("Union Exon Gene GC Content") +
  ylab("Biomart Longest Isoform GC Content") +
  ggtitle(paste0("Calc_Union_Exon_Length_and_GC.R"
                 ,"\nCompare Union Exon Gene GC Content to Biomart Longest Isoform"
                 ,"\nPearson:", round(cor(df$UnionGCcontent
                                          , df$percentage_gc_content/100
                                          , method = "pearson"
                                          , use = 'pairwise.complete.obs'), 2)
                 ,"\nSpearman:", round(cor(df$UnionGCcontent
                                           , df$percentage_gc_content/100
                                           , method = "spearman"
                                           , use = 'pairwise.complete.obs'), 2)))
ggsave("../analysis/Calc_Union_Exon_Length_and_GC_Compare_GC_To_Biomart.pdf")
################################################################################

### Calculate average length and gc for each sample

unionGenes <- data.frame(Length = ENSEMBLhg19.70UnionAnno[ ,1])
exLenDF <- merge(x = exDatDF, y = unionGenes, by.x = "X", by.y = "row.names" )
avgLength <- apply(exLenDF, 2
                 , function(counts) sum(as.numeric(counts) * exLenDF["Length"]) / 
                   sum(as.numeric(counts)))
avgLength <- tail(head(avgLength, -1), -1)
avgGCdF <- merge(x = exDatDF, y = ENSEMBLhg19.70UnionAnno, by.x = "X", by.y = "row.names" )
avgGCdF <- avgGCdF[complete.cases(avgGCdF), ]
avgGC <- apply(avgGCdF, 2
                 , function(counts) sum(as.numeric(counts) * avgGCdF["UnionGCcontent"]) / 
                   sum(as.numeric(counts)))
avgGC <- tail(head(avgGC, -2), -1)
save(avgLength, avgGC, file = "../analysis/tables/Avg_Gene_Length_and_GC.rda")
