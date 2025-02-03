# TODO:
# fragilityRow & fragilityRowNormalized
# ridge
# ridgeR2
# ridgesearchlambdadichomotomy

# setup ------------------------------------------------------------------------ 
set.seed(123L)
mat <- matrix(rnorm(1e3L, -10, 10), ncol = 10L)
ARGS = ERR <- list(xt = mat[1L:4L, ], xtp1 = mat[1L:4L + 1L, ], 0.1)
ERR$xtp <- ERR$xtp[, -1L]
RSLD <- "ridgesearchlambdadichomotomy"

# fragilityRowNormalized -------------------------------------------------------
test_that("fragilityRow & fragilityRowNormalized", {
  input <- do.call(ridge, ARGS)
  fragilityRow(input) |> expect_no_error()
  fragilityRowNormalized(input) |> expect_no_error()
})

# ridge ------------------------------------------------------------------------
test_that("ridge", {
  do.call(ridge, ERR)  |> expect_error()
  do.call(ridge, ARGS) |> expect_no_error()
})

# ridgeR2 ----------------------------------------------------------------------
test_that("ridgeR2", {
  input <- do.call(ridge, ARGS) |> list() |> c(ARGS[-3L])
  do.call(ridgeR2, input) |> expect_no_error()
})

# ridgesearchlambdadichomotomy -------------------------------------------------
test_that(RSLD, {
  c( ERR[-3L], FALSE) |> do.call(what = RSLD) |> expect_error()
  c(ARGS[-3L], TRUE)  |> do.call(what = RSLD) |> expect_no_error()
})







