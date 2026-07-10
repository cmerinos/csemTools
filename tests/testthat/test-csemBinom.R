test_that("csemBinom returns valid dataframe", {

  # Arrange
  res <- csemBinom(
    score.type = "dich",
    nitems = 16,
    ci = TRUE,
    csem.method = "Lord"
  )

  # Assert
  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) > 0)
  expect_true(ncol(res) >= 2)
})

test_that("csemBinom runs without error", {

  expect_no_error(
    csemBinom(
      score.type = "dich",
      nitems = 16,
      ci = TRUE,
      csem.method = "Lord"
    )
  )
})
