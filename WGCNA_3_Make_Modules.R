# Cluster samples and construct modules testing different parameters

# Workflow
#   1a-allen-subset-to-basal-ganglia.R
#   1b-allen-combine-probes.R
#   2a-allen-soft-thresholding-power.R
#   2b-allen-adjacency-TOM.R
#   3-allen-construct-network-modules.R
#   4-allen-compare-modules-metadata.R

print("#######################################################################")
print("Starting allen-construct-network-modules.R script...")
sessionInfo()
rm(list=ls())

library(WGCNA)
library(cluster)
library(flashClust)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

load("../data/WGCNA_2_Adjacency_TOM.rda")
load("../data/Expr_Ft_1E10S.rda")
exDatDF <- t(exDatFiltDF)

softPower <- 17

geneTree <- hclust(as.dist(dissTOM), method = "average")

# Module identification using hybrid tree cut:
# Function to construct modules
# Args: minimum module size, deep split
Make_modules <- function (minModSize, deepSplit) {
  print("Treecut arguments:")
  # print(c("minModSize"=minModSize,"cut.height"=cutHeightMergeME, "deepSplit"=ds))
  print(c(minModSize, deepSplit))
  tree = cutreeHybrid(dendro = geneTree, pamRespectsDendro = FALSE
                      , minClusterSize = minModSize
                      # , cutHeight = 0.967
                      , deepSplit = deepSplit, distM = as.matrix(dissTOM))
  print("Table of genes per module:")
  print(table(tree$labels))
  tree$labels
}

# Merge modules based on ME function
# Args: Modules colors, Cut height to merge ME
Merge_modules_ME <- function (genesModuleColor, cutHeightMergeME) {
  # Call an automatic merging function
  # merged: The merged module colors
  # Cut height of 0.25, corresponds to a correlation of 0.75, to merge ME:
  merged <- mergeCloseModules(exprData = exDatDF, colors = genesModuleColor,
                              cutHeight = cutHeightMergeME)
  labels2colors(merged$colors)
}

# Test different parameters for constructing and merging modules
# Define arguments to test for cutreeHybrid
minModSizes <- c(30, 100, 160)
deepSplits <- c(2, 4)
cutHeightMergeMEs <- c(0.1, 0.2, 0.25)
modulesColors <- NULL
moduleParameterLabels <- NULL
for (minModSize in minModSizes) {
  for (deepSplit in deepSplits) {
    # Test multiple cutreeHybrid parameters
    module <- Make_modules(minModSize, deepSplit)
    for (cutHeightMergeME in cutHeightMergeMEs) {
      # Test ME merge cut heights
      modulesColors <- cbind(modulesColors,
                             Merge_modules_ME(module, cutHeightMergeME))
      # Make label from parameters used to make each module
      moduleParameterLabels <- c(moduleParameterLabels, paste(
        "MMS=",minModSize
        , " \nDS=",deepSplit
        , " \nMEcor=",cutHeightMergeME
      ))
    }
  }
}

sizeGrWindow(25,20)
pdf("../analysis/graphs/WGCNA_3_Dendro_Module_Parameters.pdf",height=25,width=20)
plotDendroAndColors(geneTree
                    , modulesColors
                    , groupLabels=moduleParameterLabels
                    , addGuide=TRUE
                    , dendroLabels=FALSE
                    , main="Dendrogram With Different Module Cutting Parameters")
dev.off()

save(exprData, geneTree, modulesColors, moduleParameterLabels,
     file="../processed_data/allen_modules.rda")
