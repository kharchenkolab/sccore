test_that("appendSpecificityMetricsToDE AUC matches pROC::auc", {
  skip_if_not_installed("pROC")
  set.seed(1)
  nc <- 200L; ng <- 12L
  cm <- Matrix::Matrix(matrix(rpois(nc * ng, 0.6), nc, ng), sparse = TRUE)
  rownames(cm) <- paste0("c", seq_len(nc))
  colnames(cm) <- paste0("g", seq_len(ng))
  clusters <- stats::setNames(factor(sample(c("A", "B"), nc, replace = TRUE)), rownames(cm))
  de.df <- data.frame(Gene = colnames(cm), stringsAsFactors = FALSE)

  res <- appendSpecificityMetricsToDE(de.df, clusters, "A", cm, append.auc = TRUE)

  ## reference: pROC::auc per gene, fixed up-regulation direction (cases = in-cluster have higher
  ## predictor) -- this is the meaningful, presence-oriented direction the closed form computes.
  mask <- as.integer(clusters == "A")
  cb <- cm; cb@x <- as.numeric(cb@x > 0)
  ref <- vapply(colnames(cm), function(g) {
    as.numeric(pROC::auc(mask, as.integer(cb[, g]), quiet = TRUE, direction = "<", levels = c(0, 1)))
  }, numeric(1))

  expect_equal(unname(res$AUC), unname(ref), tolerance = 1e-8)
})
