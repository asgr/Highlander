test_that("Lowlander_kept ranks finite samples in the recorded direction", {
  result <- list(
    latin_mod = matrix(c(0, 10, 1, 11, 2, 12, 3, 13), ncol = 2, byrow = TRUE),
    keep = 1:4,
    output = c(4, 1, 3, 2),
    liketype = "min"
  )
  result$output[2] <- Inf

  expect_equal(
    Lowlander_kept(result),
    matrix(c(3, 13, 2, 12, 0, 10), ncol = 2, byrow = TRUE)
  )
})
