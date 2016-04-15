
# Damon Polioudakis
# 2016-04-07
# Add some RNA STAR statistics and Picard Statistics to sample metadata
# Gene Length and GC bias were calculated with Calc_Gene_Length_and_GC.R
################################################################################

rm(list=ls())
sessionInfo()

library(xlsx)

# Load data and assign variables

# Picard Sequencing Statistics
picStatsDF <- read.csv("../metadata/PicardToolsQC.csv")
dim(picStatsDF)

# RNA STAR Stats
stStatsDF <- read.table("../metadata/RNAstar_Stats.txt"
                         , sep = "\t", header = TRUE, fill = TRUE)
dim(stStatsDF)

# Metadata
metDatDF <- read.xlsx("../metadata/Cell paper - updated attributes.xlsx", 1)
dim(metDatDF)

# GC content and average length for each sample
load("../analysis/tables/Avg_Gene_Length_and_GC.rda")
avgLengthDF <- data.frame(AvgGeneLength = avgLength)
dim(avgLengthDF)
avgGCdF <- data.frame(GCcontent = avgGC)
dim(avgGCdF)
################################################################################

### Format data and add to metadata

##
# Calculate and add percent duplication data to metadata table

# Calculate percent duplication
pctDupDF <- data.frame(Pct_Duplicates =
                  picStatsDF$READ_PAIR_DUPLICATES / picStatsDF$PF_READS_ALIGNED
                  , row.names = picStatsDF$X)
# And add percent duplication data to metadata table
metDatDF <- merge(x = metDatDF, y = pctDupDF, by.x = "Cell", by.y = "row.names")

##
# Calculate and add total reads, percent uniquely mapped reads, percent mapped
# to multiple loci, and percent unmapped reads to metadata table

# Total Reads
metDatDF <- merge(x = metDatDF
                  , y = stStatsDF[ ,c("SampleID", "Number.of.input.reads")]
                  , by.x = "Cell", by.y = "SampleID")
colnames(metDatDF)[13] <- "Total_Reads"

# Uniquely Mapped - Convert % to numeric and multiply by number of reads in lane
umReadsDF <- data.frame(SampleID = stStatsDF$SampleID
    , Uniquely_Mapped = (((as.numeric(sub("%", "", stStatsDF$Uniquely.mapped.reads..))/100) * stStatsDF$Number.of.input.reads)))
metDatDF <- merge(x = metDatDF, y = umReadsDF, by.x = "Cell", by.y = "SampleID")
# Percent Uniquely Mapped
metDatDF$Pct_Uniquely_Mapped <- (metDatDF$Uniquely_Mapped / metDatDF$Total_Reads)

# Total mapped to multiple loci
metDatDF <- merge(x = metDatDF
     , y = stStatsDF[ ,c("SampleID", "Number.of.reads.mapped.to.multiple.loci")]
     , by.x = "Cell", by.y = "SampleID")
colnames(metDatDF)[16] <- "Mapped_Multiple_Loci"
# Percent mapped to multiple loci
metDatDF$Pct_Mapped_Multiple_Loci <- metDatDF$Mapped_Multiple_Loci / metDatDF$Total_Reads

# Percent unmapped reads
Calc_Unmapped <- function (data) {
  # Convert % to numeric and multiply by number of reads in lane
  df <- apply(data[ ,c("X..of.reads.unmapped..too.many.mismatches"
                                      , "X..of.reads.unmapped..too.short"
                                      , "X..of.reads.unmapped..other")]
              , 2, function(column) {as.numeric(sub("%", "", column)) / 100})
  # Multiple by number of reads in lane to convert to read number
  df <- apply(df, 2, function(column) {column * data$Number.of.input.reads})
  apply(df, 1, sum)
  }
# Calc percent unmapped
unMap <- Calc_Unmapped(stStatsDF)
unMapDF <- data.frame(SampleID = stStatsDF$SampleID, Un_Mapped = unMap)
metDatDF <- merge(x = metDatDF, y = unMapDF, by.x = "Cell", by.y = "SampleID")
# Percent unmapped reads
metDatDF$Pct_Unmapped <- metDatDF$Un_Mapped / metDatDF$Total_Reads

##
# Add average gene length to metadata

row.names(avgLengthDF) <- gsub("\\.", "_", row.names(avgLengthDF))
metDatDF <- merge(x = metDatDF, y = avgLengthDF, by.x = "Cell", by.y = "row.names")

##
# Add average GC content to metadata

row.names(avgGCdF) <- gsub("\\.", "_", row.names(avgGCdF))
metDatDF <- merge(x = metDatDF, y = avgGCdF, by.x = "Cell", by.y = "row.names")

##
# Write out to tab separated table
head(metDatDF)
write.table(metDatDF, file = paste0("../metadata/Compiled_Metadata_"
                                , format(Sys.time(), "%Y%m%d"), ".txt")
          , quote = FALSE, row.names = FALSE, sep = "\t")
