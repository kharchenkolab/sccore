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

test_that("avg_rank matches base average ranks (with ties)", {
  x <- c(3, 1, 1, 2, 5, 2, 0, 0)
  expect_equal(avg_rank(x), rank(x, ties.method = "average"))
})

test_that("sparse_matrix_column_ranks matches dense average ranks at nonzero entries", {
  set.seed(1)
  m <- matrix(rpois(10 * 6, 0.7), 10, 6)
  sm <- methods::as(Matrix::Matrix(m, sparse = TRUE), "dgCMatrix")
  rfull <- as.matrix(sparse_matrix_column_ranks(sm))   # ranks stored at nonzero, 0 elsewhere
  dense <- apply(m, 2, rank, ties.method = "average")  # full average ranks incl. zeros
  nz <- m != 0
  expect_equal(rfull[nz], dense[nz])
})

test_that("matrixDE recovers planted markers (Z, highest)", {
  set.seed(42)
  nc <- 300L; ng <- 40L
  base <- matrix(stats::rpois(nc * ng, 0.5), nc, ng)
  grp <- factor(rep(c("A", "B", "C"), length.out = nc))
  base[grp == "A", 1] <- base[grp == "A", 1] + stats::rpois(sum(grp == "A"), 5)  # g1 up in A
  base[grp == "B", 2] <- base[grp == "B", 2] + stats::rpois(sum(grp == "B"), 5)  # g2 up in B
  cm <- methods::as(Matrix::Matrix(base, sparse = TRUE), "dgCMatrix")
  rownames(cm) <- paste0("c", seq_len(nc)); colnames(cm) <- paste0("g", seq_len(ng))
  names(grp) <- rownames(cm)

  de <- matrixDE(cm, grp, z.threshold = 0)
  expect_setequal(names(de), levels(grp))
  expect_gt(de[["A"]]["g1", "Z"], 3)
  expect_true(de[["A"]]["g1", "highest"])
  expect_gt(de[["B"]]["g2", "Z"], 3)
  expect_true(de[["B"]]["g2", "highest"])

  ## specificity metrics can be appended in one call
  de2 <- matrixDE(cm, grp, z.threshold = 0, append.specificity.metrics = TRUE, append.auc = TRUE)
  expect_true(all(c("AUC", "Specificity", "Precision", "ExpressionFraction") %in% colnames(de2[["A"]])))
})
