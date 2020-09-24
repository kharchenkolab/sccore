library(sccore)
library(testthat)

test_that("setMinMax() functionality", {
	example_matrix =  matrix(rep(c(1:5), 3), 5)
	expect_equal(unique(setMinMax(example_matrix, 2, 4)[1,]), 2)
})


test_that("plapply() functionality", {
	square = function(x){ x**2 }
	expect_equal( plapply(1:10, square, n.cores=1, progress=FALSE)[[10]], 100)
})


test_that("splitVectorByNodes() functionality", {
	adjList = graphToAdjList(conosGraph)
	expect_equal(length(names(adjList)), 3)
	expect_equal(length(adjList$names), 12000)
})


test_that("graphToAdjList() functionality", {
	edge.list.fact <- igraph::as_edgelist(conosGraph) %>% sccore:::as_factor()
	edge.list <- matrix(edge.list.fact$values, ncol=2)
	n.nodes <- length(igraph::V(conosGraph))
	splitVecs = splitVectorByNodes(edge.list[,1], edge.list[,2], n.nodes)
	expect_equal(length(splitVecs), 12000)
})


test_that("graphToAdjList() functionality", {
	edge.list.fact <- igraph::as_edgelist(conosGraph) %>% sccore:::as_factor()
	edge.list <- matrix(edge.list.fact$values, ncol=2)
	n.nodes <- length(igraph::V(conosGraph))
	splitVecs = splitVectorByNodes(edge.list[,1], edge.list[,2], n.nodes)
	expect_equal(length(splitVecs), 12000)
})


test_that("graphToAdjList() functionality", {
	edge.list.fact <- igraph::as_edgelist(conosGraph) %>% sccore:::as_factor()
	edge.list <- matrix(edge.list.fact$values, ncol=2)
	n.nodes <- length(igraph::V(conosGraph))
	splitVecs = splitVectorByNodes(edge.list[,1], edge.list[,2], n.nodes)
	expect_equal(length(splitVecs), 12000)
})


test_that("embedKnnGraph() functionality", {
	adjList = graphToAdjList(conosGraph)
	commuteTimes = sccore:::get_nearest_neighbors(adjList$idx, adjList$probabilities, min_prob=1e-3, min_visited_verts=1000, n_cores=1, max_hitting_nn_num=0, max_commute_nn_num=0, min_prob_lower=1e-5, verbose=TRUE)
	knnGraph = embedKnnGraph(commuteTimes, n.neighbors=40, n.cores=1, n.epochs=1000, spread=15, min.dist=0.001, verbose=TRUE)
	expect_equal(length(knnGraph), 24000)
})


test_that("embedGraphUmap() functionality", {
	umapEmbedding = embedGraphUmap(conosGraph, verbose=TRUE, return.all=FALSE, n.cores=2)
	expect_equal(length(umapEmbedding), 24000)
})











