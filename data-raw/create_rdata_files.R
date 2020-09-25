
## Following: https://github.com/kharchenkolab/conos/blob/master/vignettes/walkthrough.md

library(conos)
library(dplyr)

panel <- readRDS(file.path(find.package('conos'),'extdata','panel.rds'))

library(pagoda2)
panel.preprocessed <- lapply(panel, basicP2proc, n.cores=4, min.cells.per.gene=0, n.odgenes=2e3, get.largevis=FALSE, make.geneknn=FALSE)

con <- Conos$new(panel.preprocessed, n.cores=4)

conosClusterList = lapply(con$samples, conos:::getCountMatrix, transposed=TRUE)
## subset the list to be smaller for the R package
conosClusterList = list(head(conosClusterList$MantonBM1_HiSeq_1), head(conosClusterList$MantonBM2_HiSeq_1))
save(conosClusterList, file="conosClusterList.rda")

con$buildGraph(k=30, k.self=5, space='PCA', ncomps=30, n.odgenes=2000, matching.method='mNN', metric='angular', score.component.variance=TRUE, verbose=TRUE)

con$findCommunities(method =igraph::walktrap.community, steps=7)

con$embedGraph(method="UMAP", min.dist=0.01, spread=15, n.cores=4, min.prob.lower=1e-3)

## conosGraph.rda
conosGraph = con$graph
save(conosGraph, file="conosGraph.rda")

umapEmbedding = conos:::embedGraphUmap(con$graph, verbose=TRUE, return.all=FALSE, n.cores=2)

## umapEmbedding.rda
save(umapEmbedding, file="umapEmbedding.rda")

cellannot <- read.table(file.path(find.package('conos'),'extdata','cellannot.txt'),header=FALSE,sep='\t')
cellAnnotations <- setNames(cellannot[,2], cellannot[,1])

## cellAnnotations.rda
save(cellAnnotations, file="cellAnnotations.rda")


