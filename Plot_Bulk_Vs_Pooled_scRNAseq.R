
# Damon Polioudakis
# 2016-04-07
# Plot bulk RNAseq of VZ and CP from Luis and Jason's ATAC versus pooled
# scRNAseq VZ from Kriegstein 2015

################################################################################

rm(list=ls())
sessionInfo()

library(ggplot2)
library(reshape2)

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

# Bulk metadata
metDatBulkDF <- read.csv("../metadata/VZCP_sampleinfo.csv", header = TRUE)

## Variables

graphCodeTitle <- "Plot_Bulk_Vs_Pooled_scRNAseq.R"
outGraph <- "../analysis/graphs/Plot_Bulk_Vs_Pooled_scRNAseq_"

################################################################################

### Format and Filter

# Select VZ samples (Kriegstein samples are all VZ)
buExDatDF <- buExDatDF[ ,metDatBulkDF$ExpCondition == "VZ"]

# Remove ERCCs
buExDatDF <- head(buExDatDF, -97)
scExDatDF <- head(scExDatDF, -97)
################################################################################

# Average counts for bulk for each transcript and counts of scRNAseq for each
# transript

# Pool scRNAseq
pScEx <- apply(scExDatDF, 1, sum)
head(pScEx)
tail(pScEx)

# Down sample to 130 cells (same number as our good cells from C196-001-002 run)
set.seed(11)
rNums <- sample(1:130, 130, replace = FALSE)
pSc130Ex <- apply(scExDatDF[ ,rNums], 1, sum)

# Mean expression for bulk for each gene
mnBuEx <- apply(buExDatDF, 1, mean)

# Median for bulk for each gene
mdBuEx <- apply(buExDatDF, 1, median)
################################################################################

## Log2 (counts + 1) Pooled scRNAseq vs bulk

# Spearman correlation
sprCorMn <- round(cor(pScEx, mnBuEx, method = "spearman"), 2)

# Format for ggplot2
ggDF <- data.frame(Pooled = log(pScEx + 1, 2), Bulk = log(mnBuEx + 1, 2))

ggplot(ggDF, aes(x = Pooled, y = Bulk)) +
  geom_point(alpha = 0.3, shape = 1) +
  theme_bw(base_size = 18) +
  ylab("Bulk: log2(Mean Counts + 1)") +
  xlab("Pooled: log2(Counts + 1)") +
  ggtitle(paste0(graphCodeTitle
                 , "\nPooled Kriegstein 2015 scRNAseq (VZ) vs Bulk RNAseq (VZ)"
                 , "\nHuman Fetal Brain"
                 , "\n Mean of counts across samples"
                 , "\nSpearman correlation: ", sprCorMn))
ggsave(paste0(outGraph, "Bulk_vs_Pooled_log2p1.pdf"))

## Down sampled to 130 cells (same number as our good cells from C196-001-002 run)
# Format for ggplot2
ggDF <- data.frame(Pooled = log(pSc130Ex + 1, 2), Bulk = log(mnBuEx + 1, 2))

ggplot(ggDF, aes(x = Pooled, y = Bulk)) +
  geom_point(alpha = 0.3, shape = 1) +
  theme_bw(base_size = 18) +
  ylab("Bulk: log2(Mean Counts + 1)") +
  xlab("Pooled: log2(Counts + 1)") +
  ggtitle(paste0(graphCodeTitle
                 , "\nPooled Kriegstein 2015 scRNAseq (VZ) vs Bulk RNAseq (VZ)"
                 , "\nDown sampled to 130 cells (same number as our good cells from C196-001-002 run)"
                 , "\nHuman Fetal Brain"
                 , "\n Mean of counts across samples"
                 , "\nSpearman correlation: ", sprCorMn))
ggsave(paste0(outGraph, "Bulk_vs_Pooled_130cells_log2p1.pdf"))
################################################################################

### Filter, read depth normalize, and pool in groups bulk and scRNAseq

## Filter for number of mapped reads > 1.5*10^6
ftScExDatDF <- scExDatDF
# # Identify Cells IDs with number of mapped reads > 1.5*10^6
# pfCells <- picStatsScDF$X[picStatsScDF$PF_READS_ALIGNED > 1.5*10^6]
# # Filter expression dataframe for cell IDs with number of mapped reads > 1.5*10^6
# ftScExDatDF <- scExDatDF[ ,colnames(scExDatDF) %in% pfCells]


## Read depth normalize

# Randomly split into 5 pooled groups of 26 cells each

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

# Mean counts for pooled scRNAseq groups
mnPdScEx <- apply(pScTpmExDatDF, 1, mean)

# Spearman correlation
sprCorMn <- round(cor(mnPdScEx, mnBuEx, method = "spearman"), 2)

# Mean counts (TPM from number mapped to Exons)
mnBuTpmExonEx <- apply(buTpmExonExDatDF, 1, mean)
mnPdScTpmExonEx <- apply(pScTpmExonExDatDF, 1, mean)


## Median counts

# Median counts for bulk RNAseq
mdBuEx <- apply(buTpmExDatDF, 1, median)

# Median counts for pooled scRNAseq groups
mdPdScEx <- apply(pScTpmExDatDF, 1, median)

# Spearman correlation
sprCorMd <- round(cor(mdPdScEx, mdBuEx, method = "spearman"), 2)
################################################################################

### Plot read depth normalized

# Boxplot read depth normalized bulk log2 (counts)
boxplot(log(data.frame(buTpmExDatDF, pScTpmExDatDF) + 1, 2), range = 0)


## Log2 mean bulk vs pooled TPM

# Format for ggplot2
ggDF <- data.frame(Pooled = log(mnPdScEx + 1, 2), Bulk = log(mnBuEx + 1, 2))

ggplot(ggDF, aes(x = Pooled, y = Bulk)) +
  geom_point(alpha = 0.5, shape = 1) +
  stat_smooth() +
  theme_bw(base_size = 18) +
  ylab("Bulk: log2(Mean TPM + 1)") +
  xlab("Pooled: log2(Mean TPM + 1)") +
  ggtitle(paste0(graphCodeTitle
                 , "\nFive Pools Kriegstein 2015 scRNAseq (VZ) vs Bulk RNAseq (VZ)"
                 , "\nHuman Fetal Brain"
                 , "\nMean of TPM across samples"
                 , "\nSpearman correlation: ", sprCorMn))
ggsave(paste0(outGraph, "Bulk_vs_Pools_TPM_log2p1.pdf"))


## Log2 median bulk vs pooled TPM

# Format for ggplot2
ggDF <- data.frame(Pooled = log(mdPdScEx + 1, 2), Bulk = log(mdBuEx + 1, 2))

ggplot(ggDF, aes(x = Pooled, y = Bulk)) +
  geom_point(alpha = 0.5, shape = 1) +
  stat_smooth() +
  theme_bw(base_size = 18) +
  ylab("Bulk: log2(Median TPM + 1)") +
  xlab("Pooled: log2(Median TPM + 1)") +
  ggtitle(paste0(graphCodeTitle
                 , "\nFive Pools Kriegstein 2015 scRNAseq (VZ) vs Bulk RNAseq (VZ)"
                 , "\nHuman Fetal Brain"
                 , "\nMedian of TPM across samples"
                 , "\nSpearman correlation: ", sprCorMd))
ggsave(paste0(outGraph, "Bulk_vs_Pools_TPM_median_log2p1.pdf"))


## MA Plot TPM (mapped reads)

# Added + 0.01 to all counts to prevent Inf values
ggDF <- data.frame(Avg = (0.5*log((mnPdScEx + 0.01) * (mnBuEx + 0.01), 2))
                   , Log2Ratio = log(((mnPdScEx + 0.01) / (mnBuEx + 0.01)), 2))

head(ggDF)
ggplot(ggDF, aes(x = Avg, y = Log2Ratio)) +
  geom_point(shape = 1, alpha = 0.5) +
  stat_smooth() +
  geom_vline(xintercept = 0, col = "red") +
  geom_hline(yintercept = -0.5, col = "red") +
  geom_hline(yintercept = 0.5, col = "red") +
  theme_bw(base_size = 14) +
  xlab("0.5*log2((Mean Pools TPM + 0.01) * (Mean Bulk TPM + 0.01))") +
  ylab("log2((Mean Pools TPM + 0.01) / (Mean Bulk TPM + 0.01))") +
  ggtitle(paste0(graphCodeTitle
                 , "\nMA Plot: Five Pools Kriegstein 2015 scRNAseq (VZ) vs Bulk RNAseq (VZ)"
                 , "\nHuman Fetal Brain"
                 , "\nMean of TPM across samples"))
ggsave(paste0(outGraph, "Bulk_vs_Pools_TPM_MAPlot_TPM_mean_p1.pdf"))


## MA Plot: Mean Average > 3.5 (~10 TPM)

# Added +1 to all counts to prevent Inf values
ggDF <- data.frame(Avg = (0.5*log((mnPdScEx + 1) * (mnBuEx + 1), 2))
                   , Log2Ratio = log(((mnPdScEx + 1) / (mnBuEx + 1)), 2))
ggDF <- ggDF[ggDF$Avg > 3.5, ]

head(ggDF)
ggplot(ggDF, aes(x = Avg, y = Log2Ratio)) +
  geom_point(shape = 1, alpha = 0.5) +
  stat_smooth() +
  theme_bw(base_size = 14) +
  xlab("0.5*log2((Mean Pools TPM + 1) * (Mean Bulk TPM + 1))") +
  ylab("log2((Mean Pools TPM + 1) / (Mean Bulk TPM + 1))") +
  ggtitle(paste0(graphCodeTitle
                 , "\nMA Plot: Five Pools Kriegstein 2015 scRNAseq (VZ) vs Bulk RNAseq (VZ)"
                 , "\nHuman Fetal Brain"
                 , "\nMean of TPM across samples"
                 , "\nFiltered for Mean Average > 3.5 (~10 TPM)"))
ggsave(paste0(outGraph, "Bulk_vs_Pools_TPM_MAPlot_AG3.5_TPM_mean_p1.pdf"))


## MA Plot: TPM (Mapped to Exons)

# Added + 0.01 to all counts to prevent Inf values
ggDF <- data.frame(Avg = (0.5*log((mnPdScTpmExonEx + 0.01) * (mnBuTpmExonEx + 0.01), 2))
                   , Log2Ratio = log(((mnPdScTpmExonEx + 0.01) / (mnBuTpmExonEx + 0.01)), 2))

head(ggDF)
ggplot(ggDF, aes(x = Avg, y = Log2Ratio)) +
  geom_point(shape = 1, alpha = 0.5) +
  stat_smooth() +
  geom_vline(xintercept = 0, col = "red") +
  geom_hline(yintercept = -0.5, col = "red") +
  geom_hline(yintercept = 0.5, col = "red") +
  theme_bw(base_size = 14) +
  xlab("0.5*log2((Mean Pools TPM + 0.01) * (Mean Bulk TPM + 0.01))") +
  ylab("log2((Mean Pools TPM + 0.01) / (Mean Bulk TPM + 0.01))") +
  ggtitle(paste0(graphCodeTitle
                 , "\nMA Plot: Five Pools Kriegstein 2015 scRNAseq (VZ) vs Bulk RNAseq (VZ)"
                 , "\nHuman Fetal Brain"
                 , "\nMean of TPM (read depth normalized by number mapped to exons) across samples"))
ggsave(paste0(outGraph, "Bulk_vs_Pools_TPM_MAPlot_TPMexon_mean_p1.pdf"))

# stat_smooth(method = lm)
# + coord_cartesian(xlim = c(0, 500000), ylim = c(0, 2500))
################################################################################

### CoV versus TPM

# Artifactual looking points at top of ScCov are most likely due to 4 of the
# pooled groups having 0 counts and 1 group having 1 count, not sure what the
# 1 normalizes to after TPM, but that math should carry through

## TPM
rDep <- (apply(ftScExDatDF, 2, sum) / 10^6)
scRdnExDatDF <- ftScExDatDF / rDep


## Mean TPM

# Mean TPM for bulk RNAseq
mnBuEx <- apply(buTpmExDatDF, 1, mean)
lg2mnBuEx <- log(mnBuEx, 2)

# Mean TPM for pooled scRNAseq groups
mnPdScEx <- apply(pScTpmExDatDF, 1, mean)
lg2mnPdScEx <- log(mnPdScEx, 2)

# Mean TPM for scRNAseq
mnScEx <- apply(scRdnExDatDF, 1, mean)
lg2mnScEx <- log(mnScEx, 2)


## CoV

Coef_Of_Var <- function(x) {
  sd(x) / mean(x)
}

buCov <- apply(buTpmExDatDF, 1, function(tpm) Coef_Of_Var(tpm))
pdScCov <- apply(pScTpmExDatDF, 1, function(tpm) Coef_Of_Var(tpm))
scCov <- apply(scRdnExDatDF, 1, function(tpm) Coef_Of_Var(tpm))


## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 18)))
theme_update(plot.title = element_text(size = 14))


## Plot CoV as panels
ggDF <- data.frame(Bulk = buCov, Pooled = pdScCov, scRNAseq = scCov)
ggDF <- melt(ggDF)
ggDF$Mean_lg2TPM <- c(lg2mnBuEx, lg2mnPdScEx, lg2mnScEx)
nMisVals <- summary(sapply((split(ggDF$value, ggDF$variable)), is.na))
nMisVals
ggplot(ggDF, aes(x = Mean_lg2TPM, y = value)) +
  geom_point(shape = 1, alpha = 0.25, size = 0.5) +
  facet_wrap(~variable, scales = "free") +
  xlab("Log2(Mean TPM)") +
  ylab("Coefficient of Variation of TPM") +
  ggtitle(paste0(graphCodeTitle
                 , "\nlog2 TPM vs CoV: RNAseq - Human Fetal Brain VZ"
                 , "\nMean of TPM across samples (TPM: normalized by number mapped reads)"
                 , "\n"))
ggsave(paste0(outGraph, "CoV_TPM.pdf"))


## Plot Bulk CoV vs TPM
ggDF <- data.frame(CoV = buCov, Mean_lg2TPM = lg2mnBuEx)
nMisVals <- table(is.na(ggDF$CoV))
ggplot(ggDF, aes(x = Mean_lg2TPM, y = CoV)) +
  geom_point(shape = 1, alpha = 0.25) +
  xlab("Log2(Mean TPM)") +
  ylab("Coefficient of Variation of TPM") +
  ggtitle(paste0(graphCodeTitle
                 , "\nlog2 TPM vs CoV: RNAseq - Human Fetal Brain VZ"
                 , "\nMean of TPM across samples (TPM: normalized by number mapped reads)"
                 , "\nNAs removed after CoV calculation: "
                 , nMisVals[2], " out of ", sum(nMisVals[1:2]), " genes"
                 , "\n"))
ggsave(paste0(outGraph, "CoV_Bulk_TPM.pdf"))


## Plot pooled scRNAseq CoV vs TPM
ggDF <- data.frame(CoV = pdScCov, Mean_lg2TPM = lg2mnPdScEx)
nMisVals <- table(is.na(ggDF$CoV))
ggplot(ggDF, aes(x = Mean_lg2TPM, y = CoV)) +
  geom_point(shape = 1, alpha = 0.25) +
  xlab("Log2(Mean TPM)") +
  ylab("Coefficient of Variation of TPM") +
  ggtitle(paste0(graphCodeTitle
                 , "\nlog2 TPM vs CoV: Pooled Kriegstein 2015 scRNAseq - Human Fetal Brain VZ"
                 , "\nMean of TPM across samples (TPM: normalized by number mapped reads)"
                 , "\nNAs removed after CoV calculation: "
                 , nMisVals[2], " out of ", sum(nMisVals[1:2]), " genes"
                 , "\n"))
ggsave(paste0(outGraph, "CoV_pooledSc_TPM.pdf"))


## Plot scRNAseq CoV vs TPM
ggDF <- data.frame(CoV = scCov, Mean_lg2TPM = lg2mnScEx)
nMisVals <- table(is.na(ggDF$CoV))
ggplot(ggDF, aes(x = Mean_lg2TPM, y = CoV)) +
  geom_point(shape = 1, alpha = 0.25) +
  xlab("Log2(Mean TPM)") +
  ylab("Coefficient of Variation of TPM") +
  ggtitle(paste0(graphCodeTitle
                 , "\nlog2 TPM vs CoV: Kriegstein 2015 scRNAseq - Human Fetal Brain VZ"
                 , "\nMean of TPM across samples (TPM: normalized by number mapped reads)"
                 , "\nNAs removed after CoV calculation: "
                 , nMisVals[2], " out of ", sum(nMisVals[1:2]), " genes"
                 , "\n"))
ggsave(paste0(outGraph, "CoV_scRNAseq_TPM.pdf"))


### Median absolute deviation / median rather than coefficient of variation

# The steps to find the MAD include:
# 1. find the mean (average)
# 2. find the difference between each data value and the mean
# 3. take the absolute value of each difference
# 4. find the mean (average) of these differences

## MAD / median
MADmedian <- function(x) {
  tmp <- sapply(x, function(tpm) abs(tpm - median(x)))
  median(tmp) / median(x)
}
buMADm <- apply(buTpmExDatDF, 1, function(tpm) MADmedian(tpm))
pdScMADm <- apply(pScTpmExDatDF, 1, function(tpm) MADmedian(tpm))
scMADm <- apply(scRdnExDatDF, 1, function(tpm) MADmedian(tpm))
# High number of 0 values results from MAD frequently having same value as median


## Median counts

# Median counts for bulk RNAseq
mdBuEx <- apply(buTpmExDatDF, 1, median)
lg2mdBuEx <- log(mdBuEx, 2)

# Median counts for pooled scRNAseq groups
mdPdScEx <- apply(pScTpmExDatDF, 1, median)
lg2mdPdScEx <- log(mdPdScEx, 2)

# Median counts for scRNAseq
mdScEx <- apply(scRdnExDatDF, 1, median)
lg2mdScEx <- log(mdScEx, 2)


## Plot MAD/median as panels
ggDF <- data.frame(Bulk = buMADm, Pooled = pdScMADm, scRNAseq = scMADm)
ggDF <- melt(ggDF)
ggDF$Median_lg2TPM <- c(lg2mdBuEx, lg2mdPdScEx, lg2mdScEx)
nMisVals <- summary(sapply((split(ggDF$value, ggDF$variable)), is.na))
nMisVals
ggplot(ggDF, aes(x = Median_lg2TPM, y = value)) +
  geom_point(shape = 1, alpha = 0.25, size = 0.5) +
  facet_wrap(~variable, scales = "free") +
  xlab("Log2(Median TPM)") +
  ylab("MAD / median of TPM") +
  ggtitle(paste0(graphCodeTitle
                 , "\nlog2 TPM vs MAD / median: RNAseq - Human Fetal Brain VZ"
                 , "\nMedian of TPM across samples (TPM: normalized by number mapped reads)"
                 , "\n"))
ggsave(paste0(outGraph, "MADmedian_TPM.pdf"))


## Plot Bulk MAD/median vs TPM
ggDF <- data.frame(MADmedian = buMADm, Median_lg2TPM = lg2mdBuEx)
nMisVals <- table(is.na(ggDF$MADmedian))
ggplot(ggDF, aes(x = Median_lg2TPM, y = MADmedian)) +
  geom_point(shape = 1, alpha = 0.25) +
  xlab("Log2(Median TPM)") +
  ylab("MAD / median of TPM") +
  ggtitle(paste0(graphCodeTitle
                 , "\nlog2 TPM vs MAD / median: Bulk RNAseq - Human Fetal Brain VZ"
                 , "\nMedian of TPM across samples (TPM: normalized by number mapped reads)"
                 , "\nNAs removed after MAD / median calculation: "
                 , nMisVals[2], " out of ", sum(nMisVals[1:2]), " genes"
                 , "\n"))
ggsave(paste0(outGraph, "MADmedian_Bulk_TPM.pdf"))


## Plot pooled scRNAseq MAD/median vs TPM
ggDF <- data.frame(MADmedian = pdScMADm, Median_lg2TPM = lg2mdPdScEx)
nMisVals <- table(is.na(ggDF$MADmedian))
ggplot(ggDF, aes(x = Median_lg2TPM, y = MADmedian)) +
  geom_point(shape = 1, alpha = 0.25) +
  xlab("Log2(Median TPM)") +
  ylab("MAD / median of TPM") +
  ggtitle(paste0(graphCodeTitle
                 , "\nlog2 TPM vs MAD / median: Pooled scRNAseq - Kriegstein 2015 Human Fetal Brain VZ"
                 , "\nMedian of TPM across samples (TPM: normalized by number mapped reads)"
                 , "\nNAs removed after MAD / median calculation: "
                 , nMisVals[2], " out of ", sum(nMisVals[1:2]), " genes"
                 , "\n"))
ggsave(paste0(outGraph, "MADmedian_pooledSc_TPM.pdf"))


## Plot scRNAseq MAD/median vs TPM
ggDF <- data.frame(MADmedian = scMADm, Median_lg2TPM = lg2mdScEx)
nMisVals <- table(is.na(ggDF$MADmedian))
ggplot(ggDF, aes(x = Median_lg2TPM, y = MADmedian)) +
  geom_point(shape = 1, alpha = 0.25) +
  xlab("Log2(Median TPM)") +
  ylab("MAD / median of TPM") +
  ggtitle(paste0(graphCodeTitle
                 , "\nlog2 TPM vs MAD / median: scRNAseq - Kriegstein 2015 Human Fetal Brain VZ"
                 , "\nMedian of TPM across samples (TPM: normalized by number mapped reads)"
                 , "\nNAs removed after MAD / median calculation: "
                 , nMisVals[2], " out of ", sum(nMisVals[1:2]), " genes"
                 , "\n"))
ggsave(paste0(outGraph, "MADmedian_scRNAseq_TPM.pdf"))