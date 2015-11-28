# Cluster samples and construct modules

# Workflow
#   1a-allen-subset-to-basal-ganglia.R
#   1b-allen-combine-probes.R
#   2a-allen-soft-thresholding-power.R
#   2b-allen-adjacency-TOM.R
#   3-allen-construct-network-modules.R
#   4-allen-compare-modules-metadata.R

print("#######################################################################")
print("Starting WGCNA_2_Calc_Adj_TOM.R script...")
sessionInfo()
rm(list=ls())

library(WGCNA)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
# disableWGCNAThreads() 

load("../data/Expr_Ft_1E10S.rda")
exDatDF <- t(exDatFiltDF)

softPower <- 17
# Biweight midcorrelation is considered to be a good alternative to Pearson
# correlation since it is more robust to outliers.

print("Starting adjacency calculation...")
adjacency <- adjacency(exDatDF, power= softPower, corFnc= "bicor")
print("Finished adjacency calculation...")

print("Starting TOM calculation...")
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM
print("Finished TOM calculation...")

save(adjacency, TOM, dissTOM
     , file = paste("../data/WGCNA_2_Adjacency_TOM_SP", softPower, ".rda"
                    , sep = ""))
print("Finished TOM calculation...")

print("End of WGCNA_2_Calc_Adj_TOM.R script...")