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

test_that("dotPlotData aggregates sparse expression by group", {
  cm <- Matrix::Matrix(
    c(
      1, 3, 0, 0,
      0, 2, 4, 0
    ),
    nrow = 4,
    ncol = 2,
    sparse = TRUE,
    dimnames = list(paste0("c", 1:4), c("g1", "g2"))
  )
  groups <- factor(c("A", "A", "B", "B"), levels = c("A", "B"))
  names(groups) <- rownames(cm)

  dp <- dotPlotData(c("g1", "g2"), cm, groups)

  expect_equal(dp$pct.exp[dp$cluster == "A" & dp$gene == "g1"], 100)
  expect_equal(dp$pct.exp[dp$cluster == "B" & dp$gene == "g1"], 0)
  expect_equal(dp$pct.exp[dp$cluster == "A" & dp$gene == "g2"], 50)
  expect_equal(dp$avg.exp[dp$cluster == "A" & dp$gene == "g1"], 2)
  expect_true(is.nan(dp$avg.exp[dp$cluster == "B" & dp$gene == "g1"]))
})

test_that("dotPlot returns a ggplot object from sparse input", {
  cm <- Matrix::Matrix(
    c(
      1, 3, 0, 0,
      0, 2, 4, 0
    ),
    nrow = 4,
    ncol = 2,
    sparse = TRUE,
    dimnames = list(paste0("c", 1:4), c("g1", "g2"))
  )
  groups <- factor(c("A", "A", "B", "B"), levels = c("A", "B"))
  names(groups) <- rownames(cm)

  expect_s3_class(dotPlot(c("g1", "g2"), cm, groups), "ggplot")
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









