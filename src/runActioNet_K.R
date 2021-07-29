args <- commandArgs(TRUE)
suppressMessages(require(ACTIONet))
suppressMessages(require(SingleCellExperiment))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))

setwd(args[1])

ace <- readRDS('ACTIONet-model')
clusters <- Leiden.clustering(ace)

write.table(as.matrix(colMaps(ace)$H_unified), file="archetypal-explicit-function.txt")	
write.table(clusters, file="assigned-subpopulation.txt")

pdf(sprintf("subpopulation-visualization.pdf"))
plot.ACTIONet(ace, clusters, node.size = 2, transparency.attr = ace$node_centrality, trans.fact = 3)
dev.off()
