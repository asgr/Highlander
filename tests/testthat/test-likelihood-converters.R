test_that("CMA converters apply optimization direction and prior consistently", {
  scalar_like <- function(parm, Data) parm[1]^2
  prior <- function(parm, Data) -0.5

  expect_equal(
    Highlander:::.convert_CMA2CMA(
      parm = 2, Data = list(), likefunc = scalar_like,
      liketype = "max", prior = prior
    ),
    -3.5
  )

  converted <- Highlander:::.convert_CMA2LD(
    parm = 2,
    Data = list(applyconstraints = FALSE, applyintervals = FALSE),
    likefunc = scalar_like,
    liketype = "min",
    prior = prior
  )
  expect_equal(converted$LP, -4.5)
  expect_equal(converted$Dev, 8)
  expect_equal(converted$parm, 2)
})

test_that("LD likelihoods are converted to the CMA objective", {
  ld_like <- function(parm, Data) list(LP = -4)

  expect_equal(
    Highlander:::.convert_LD2CMA(
      parm = 1, Data = list(), likefunc = ld_like,
      liketype = "max", prior = function(parm, Data) -0.5
    ),
    4.5
  )
})
