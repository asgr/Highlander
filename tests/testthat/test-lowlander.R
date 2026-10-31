test_that("Lowlander maps design points to bounds and selects the best sample", {
  latin <- cbind(seq(0, 1, length.out = 6), seq(0, 1, length.out = 6))
  colnames(latin) <- c("x", "y")
  result <- Lowlander(
    lower = c(0, 0),
    upper = c(1, 1),
    Data = list(),
    likefunc = function(parm, Data) sum((parm - c(0.2, 0.2))^2),
    Nsamp = nrow(latin),
    pcut = 0.75,
    latin = latin,
    parm.names = c("x", "y")
  )

  expect_equal(result$latin_mod, latin)
  expect_equal(result$best, c(x = 0.2, y = 0.2))
  expect_equal(result$keep, 1:4)
  expect_equal(result$lower, c(x = 0, y = 0))
  expect_equal(result$upper, c(x = 0.6, y = 0.6))
  expect_equal(result$frac_keep, 4 / 6)
  expect_equal(result$bad, integer())
})

test_that("Lowlander supports maximization and excludes non-finite samples", {
  latin <- cbind(seq(0, 1, length.out = 6), seq(0, 1, length.out = 6))
  result <- Lowlander(
    lower = c(0, 0),
    upper = c(1, 1),
    Data = list(),
    likefunc = function(parm, Data) {
      if (all(parm == 0.4)) return(NA_real_)
      sum(parm)
    },
    Nsamp = nrow(latin),
    pcut = 0.8,
    liketype = "max",
    latin = latin
  )

  expect_equal(result$liketype, "max")
  expect_equal(result$bad, 3L)
  expect_equal(result$best, c(1, 1))
  expect_true(all(result$keep %in% setdiff(seq_len(nrow(latin)), result$bad)))
})

test_that("Lowlander validates its design and scalar likelihood output", {
  latin <- matrix(seq(0, 1, length.out = 8), ncol = 2)
  args <- list(
    lower = c(0, 0),
    upper = c(1, 1),
    Data = list(),
    likefunc = function(parm, Data) sum(parm),
    Nsamp = nrow(latin),
    pcut = 0.5,
    latin = latin
  )

  expect_error(
    do.call(Lowlander, utils::modifyList(args, list(upper = 1))),
    "'lower' and 'upper'"
  )
  expect_error(
    do.call(Lowlander, utils::modifyList(args, list(pcut = 1))),
    "'pcut'"
  )
  expect_error(
    do.call(Lowlander, utils::modifyList(
      args, list(latin = matrix(2, nrow = nrow(latin), ncol = 2))
    )),
    "values in 'latin'"
  )
  expect_error(
    do.call(Lowlander, utils::modifyList(
      args, list(likefunc = function(parm, Data) list(LP = 1, Dev = 2))
    )),
    "single numeric value"
  )
  expect_error(
    do.call(Lowlander, utils::modifyList(
      args, list(likefunc = function(parm, Data) "not numeric")
    )),
    "must return numeric values"
  )
})
