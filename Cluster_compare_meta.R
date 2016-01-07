print("#######################################################################")
print("Starting allen-cluster-compare-meta-data.R script...")
sessionInfo()
rm(list=ls())

library(cluster)
library(flashClust)


# Cluster samples and compare to meta data

load("../data/WGCNA_3_Modules.rda")

metDatDF <- read.csv("../data/mmc2.csv", header = TRUE)

expr.data <- exDatDF

metDatDF$AlignmentRate..Pairs <- as.numeric(
  substr(as.character(metDatDF$AlignmentRate..Pairs), 1, 4))
trait.data <- metDatDF

# Cluster samples
sampleTree2 = hclust(dist(expr.data), method = "average")

# Convert traits to a color representation: for numeric traits white means low, 
# red means high, grey means missing entry
traitColors= data.frame(labels2colors(trait.data[,2:3])
                        , numbers2colors(trait.data[,4:7]))
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(trait.data)[-1],
                    main = "Sample dendrogram and trait heatmap")

print("Ending allen-cluster-compare-meta-data.R script...")