test_that("curvesCSEMcomp returns expected structure", {

  set.seed(123)

  score <- 1:20

  csem1 <- 1 + 0.01 * score + rnorm(20, 0, 0.02)
  csem2 <- 1.2 + 0.008 * score + rnorm(20, 0, 0.02)

  res <- curvesCSEMcomp(
    score = score,
    csem.m1 = csem1,
    csem.m2 = csem2
  )

  expect_s3_class(res, "csem_compare")

  expect_true("csem.compare" %in% names(res))
  expect_true("var.compare" %in% names(res))
  expect_true("global.effect" %in% names(res))
  expect_true("regression" %in% names(res))

})
