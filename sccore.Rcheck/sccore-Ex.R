pkgname <- "sccore"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "sccore-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('sccore')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("adjacentVertices")
### * adjacentVertices

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: adjacentVertices
### Title: List of adjacent vertices from igraph object
### Aliases: adjacentVertices

### ** Examples

## Not run: 
##D edges <- igraph::as_edgelist(conosGraph)
##D adjacentVertices(edges)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("adjacentVertices", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("adjacent_vertex_weights")
### * adjacent_vertex_weights

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: adjacent_vertex_weights
### Title: List of adjacent vertex weights from igraph object
### Aliases: adjacent_vertex_weights

### ** Examples

## Not run: 
##D edges <- igraph::as_edgelist(conosGraph)
##D edge.weights <- igraph::edge.attributes(conosGraph)$weight
##D adjacent_vertex_weights(edges, edge.weights)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("adjacent_vertex_weights", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("collapseGraphSum")
### * collapseGraphSum

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: collapseGraphSum
### Title: Collapse Graph By Sum
### Aliases: collapseGraphSum

### ** Examples

## No test: 
collapsed = collapseGraphPaga(conosGraph, igraph::V(conosGraph), linearize=TRUE, winsorize=FALSE)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("collapseGraphSum", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("embeddingPlot")
### * embeddingPlot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: embeddingPlot
### Title: Plot embedding with provided labels / colors using ggplot2
### Aliases: embeddingPlot

### ** Examples

library(sccore)
embeddingPlot(umapEmbedding, show.ticks=TRUE, show.labels=TRUE, title="UMAP embedding")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("embeddingPlot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("extendMatrix")
### * extendMatrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: extendMatrix
### Title: Extend matrix to include new columns in matrix
### Aliases: extendMatrix

### ** Examples

library(dplyr)
geneUnion <- lapply(conosClusterList, colnames) %>% Reduce(union, .)
extendMatrix(conosClusterList[[1]], col.names=geneUnion)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("extendMatrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fac2col")
### * fac2col

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fac2col
### Title: Utility function to translate a factor into colors
### Aliases: fac2col

### ** Examples

genes = factor(c("BRAF", "NPC1", "PAX3", "BRCA2", "FMR1"))
fac2col(genes)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fac2col", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getClusterGraph")
### * getClusterGraph

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getClusterGraph
### Title: Collapse vertices belonging to each cluster in a graph
### Aliases: getClusterGraph

### ** Examples

## No test: 
cluster.graph = getClusterGraph(conosGraph, igraph::V(conosGraph))
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getClusterGraph", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("graphToAdjList")
### * graphToAdjList

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: graphToAdjList
### Title: Convert igraph graph into an adjacency list
### Aliases: graphToAdjList

### ** Examples

library(dplyr)
edge.list.fact <- igraph::as_edgelist(conosGraph) %>% as_factor()
edge.list <- matrix(edge.list.fact$values, ncol=2)
n.nodes <- length(igraph::V(conosGraph))
splitVectorByNodes(edge.list[,1], edge.list[,2], n.nodes)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("graphToAdjList", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("jsDist")
### * jsDist

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: jsDist
### Title: Jensen<e2><80><93>Shannon distance metric (i.e. the square root
###   of the Jensen<e2><80><93>Shannon divergence) between the columns of a
###   dense matrix m
### Aliases: jsDist

### ** Examples

ex = matrix(1:9, nrow = 3, ncol = 3)
jsDist(ex)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("jsDist", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mergeCountMatrices")
### * mergeCountMatrices

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mergeCountMatrices
### Title: Merge list of count matrices into a common matrix, entering 0s
###   for the missing entries
### Aliases: mergeCountMatrices

### ** Examples

mergeCountMatrices(conosClusterList, n.cores=1)
## 12 x 67388 sparse Matrix of class "dgCMatrix"




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mergeCountMatrices", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plapply")
### * plapply

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plapply
### Title: Parallel, optionally verbose lapply. See ?parallel::mclapply for
###   more info.
### Aliases: plapply

### ** Examples

square = function(x){ x**2 }
plapply(1:10, square, n.cores=1, progress=TRUE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plapply", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("propagateLabels")
### * propagateLabels

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: propagateLabels
### Title: Estimate labeling distribution for each vertex, based on
###   provided labels.
### Aliases: propagateLabels

### ** Examples

propagateLabels(conosGraph, labels=cellAnnotations)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("propagateLabels", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("propagateLabelsDiffusion")
### * propagateLabelsDiffusion

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: propagateLabelsDiffusion
### Title: Estimate labeling distribution for each vertex, based on
###   provided labels using a Random Walk on graph
### Aliases: propagateLabelsDiffusion

### ** Examples

propagateLabelsDiffusion(conosGraph, labels=cellAnnotations)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("propagateLabelsDiffusion", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("propagateLabelsSolver")
### * propagateLabelsSolver

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: propagateLabelsSolver
### Title: Propagate labels using Zhu, Ghahramani, Lafferty (2003)
###   algorithm, "Semi-Supervised Learning Using Gaussian Fields and
###   Harmonic Functions" <http://mlg.eng.cam.ac.uk/zoubin/papers/zgl.pdf>
### Aliases: propagateLabelsSolver

### ** Examples

## No test: 
propagateLabelsSolver(conosGraph, labels=cellAnnotations)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("propagateLabelsSolver", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("setMinMax")
### * setMinMax

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: setMinMax
### Title: Set range for values in object. Changes values outside of range
###   to min or max. Adapted from Seurat::MinMax
### Aliases: setMinMax

### ** Examples

example_matrix =  matrix(rep(c(1:5), 3), 5)
setMinMax(example_matrix, 2, 4)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("setMinMax", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sn")
### * sn

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sn
### Title: Set names equal to values, a stats::setNames wrapper function
### Aliases: sn

### ** Examples

vec = c(1, 2, 3, 4)
sn(vec)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sn", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("splitVectorByNodes")
### * splitVectorByNodes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: splitVectorByNodes
### Title: splitVectorByNodes
### Aliases: splitVectorByNodes

### ** Examples

adjList = graphToAdjList(conosGraph)
print(names(adjList))
## [1] "idx" "probabilities" "names" 
length(adjList$names)
## [1] 12000




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("splitVectorByNodes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
