library(sccore)
library(testthat)

test_that("basic setMinMax() functionality", {
	example_matrix =  matrix(rep(c(1:5), 3), 5)
	expect_equal(unique(setMinMax(example_matrix, 2, 4)[1,]), 2)
})


test_that("basic plapply() functionality", {
	square = function(x){ x**2 }
	expect_equal( plapply(1:10, square, n.cores=1, progress=FALSE)[[10]], 100)
})

