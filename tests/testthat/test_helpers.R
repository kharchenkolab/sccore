library(sccore)
library(dplyr)
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
	expect_equal(length(adjList$names), 100)
})


test_that("graphToAdjList() functionality", {
	edge.list.fact <- igraph::as_edgelist(conosGraph) %>% as_factor()
	edge.list <- matrix(edge.list.fact$values, ncol=2)
	n.nodes <- length(igraph::V(conosGraph))
	splitVecs = splitVectorByNodes(edge.list[,1], edge.list[,2], n.nodes)
	expect_equal(length(splitVecs), 100)
})


test_that("graphToAdjList() functionality", {
	edge.list.fact <- igraph::as_edgelist(conosGraph) %>% as_factor()
	edge.list <- matrix(edge.list.fact$values, ncol=2)
	n.nodes <- length(igraph::V(conosGraph))
	splitVecs = splitVectorByNodes(edge.list[,1], edge.list[,2], n.nodes)
	expect_equal(length(splitVecs), 100)
})


test_that("graphToAdjList() functionality", {
	edge.list.fact <- igraph::as_edgelist(conosGraph) %>% as_factor()
	edge.list <- matrix(edge.list.fact$values, ncol=2)
	n.nodes <- length(igraph::V(conosGraph))
	splitVecs = splitVectorByNodes(edge.list[,1], edge.list[,2], n.nodes)
	expect_equal(length(splitVecs), 100)
})










