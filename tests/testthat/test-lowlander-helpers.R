test_that("Lowlander_kept filters and ranks samples in the recorded direction", {
  result <- lowlander_fixture()
  result$keep <- c(1, 2, 3, 4)
  result$output[2] <- Inf

  kept <- Lowlander_kept(result, parm.names = c("first", "second"))

  expect_equal(kept, matrix(c(2, 12, 3, 13, 0, 10), ncol = 2, byrow = TRUE,
                           dimnames = list(NULL, c("first", "second"))))

  result$output <- c(4, 1, 3, 2)
  result$liketype <- "max"
  expect_equal(Lowlander_kept(result)[1, ], c(0, 10))
  expect_null(Lowlander_kept(NULL))
})

test_that("Lowlander_Z returns spread, bounded starting points", {
  result <- lowlander_fixture()

  z <- Lowlander_Z(result, N = 2, parm.names = c("first", "second"))
  expect_equal(dim(z), c(2, 2))
  expect_equal(colnames(z), c("first", "second"))
  expect_equal(z, matrix(c(1, 11, 2, 12), ncol = 2, byrow = TRUE,
                         dimnames = list(NULL, c("first", "second"))))

  set.seed(10)
  larger_z <- Lowlander_Z(result, N = 6)
  expect_equal(dim(larger_z), c(6, 2))
  expect_true(all(larger_z[, 1] >= 0 & larger_z[, 1] <= 3))
  expect_true(all(larger_z[, 2] >= 10 & larger_z[, 2] <= 13))
  expect_equal(anyDuplicated(as.data.frame(larger_z)), 0L)
})

test_that("Lowlander helpers handle unusable inputs and name the best point", {
  expect_warning(
    expect_null(Lowlander_kept(list(keep = 1))),
    "missing latin_mod/keep/output"
  )

  invalid_keep <- lowlander_fixture()
  invalid_keep$keep <- 5L
  expect_warning(
    expect_null(Lowlander_kept(invalid_keep)),
    "no valid kept samples"
  )

  nonfinite <- lowlander_fixture()
  nonfinite$output[nonfinite$keep] <- Inf
  expect_warning(
    expect_null(Lowlander_kept(nonfinite)),
    "All kept Lowlander samples have non-finite likelihood"
  )

  best <- Lowlander_best(lowlander_fixture(), parm.names = c("first", "second"))
  expect_equal(best, c(first = 1, second = 11))
  expect_null(Lowlander_best(NULL))
  expect_warning(expect_null(Lowlander_best(list())), "no $best")
})
