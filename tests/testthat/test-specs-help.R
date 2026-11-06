test_that("Specs_help provides algorithm defaults and seeds ensemble samplers", {
  data <- list(parm.names = c("first", "second"))
  lowlander <- lowlander_fixture()

  expect_null(Specs_help("CHARM", data))

  aies <- Specs_help("AIES", data, Lowlander = lowlander, Nc = 5, CPUs = 2)
  expect_equal(aies$Nc, 6L)
  expect_equal(dim(aies$Z), c(6, 2))
  expect_true(all(aies$Z[, 1] >= lowlander$lower[1] &
                  aies$Z[, 1] <= lowlander$upper[1]))

  demc <- Specs_help("DEMC", data, Lowlander = lowlander,
                     Iterations = 12, Thinning = 2)
  expect_equal(dim(demc$Z), c(7, 2, 4))
  expect_true(all(demc$Z[, 1, ] >= lowlander$lower[1] &
                  demc$Z[, 1, ] <= lowlander$upper[1]))
})

test_that("Specs_help scales Griddy-Gibbs grids to plausible bounds", {
  data <- list(parm.names = c("first", "second"))
  specs <- Specs_help("AGG", data, Lowlander = lowlander_fixture())

  expect_length(specs$Grid, 2)
  expect_equal(specs$Grid[[1]][1], -0.75)
  expect_equal(specs$Grid[[1]][length(specs$Grid[[1]])], 0.75)
  expect_null(specs$dparm)
})
