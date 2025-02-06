# TODO:
# fragilityRow & fragilityRowNormalized
# ridge
# ridgeR2
# ridgesearchlambdadichomotomy

# setup ------------------------------------------------------------------------
set.seed(123L)
mat <- matrix(rnorm(1e3L, -10, 10), ncol = 10L)
ARGS = ERR <- list(xt = mat[1L:4L, ], xtp1 = mat[1L:4L + 1L, ], 0.1)
ERR$xtp1 <- ERR$xtp1[, -1L]
RSLD <- "ridgesearchlambdadichomotomy"

# fragilityRowNormalized -------------------------------------------------------
test_that("Fragilities", {
  input <- do.call(ridge, ARGS)
  fragilityRow(input) |> expect_no_error()
  fragilityRowNormalized(input) |> expect_no_error()
  FragilityROW(input) |> expect_no_error()
  FragilityROW(input, normalize = FALSE) |> expect_no_error()
  # capture.output({
  #   fragObj <- GetAdjFrag(matrix(rnorm(1000), 100), 50, 25, 1e-4)
  # }) |> invisible()
  # fStats <- fragStat(fragObj, 1:5) |> expect_no_error()
  # fStats@qmatrix |> expect_no_error()
  # fStats@qmatrix <- NULL |> expect_no_error()
  # fragObj@frag |> as.data.frame() |> fragStat(1:5) |> expect_error()
  # frag_stat(fragObj, 1:5) |> expect_no_error()
})

# ridge ------------------------------------------------------------------------
test_that("ridge/RIDGE", {
  do.call(ridge, ERR)  |> expect_error()
  do.call(ridge, ARGS) |> expect_no_error()
  do.call(RIDGE, unname(ARGS)) |> expect_no_error()
})

# ridgeR2 ----------------------------------------------------------------------
test_that("RidgeRSQ/ridgeR2", {
  input <- ARGS[-3L] |> c(list(do.call(ridge, ARGS))) |> unname()
  do.call(ridgeR2, input) |> expect_no_error()
  do.call(RidgeRSQ, input) |> expect_no_error()
})

# ridgesearchlambdadichomotomy -------------------------------------------------
test_that(RSLD, {
  c( ERR[-3L], FALSE) |> do.call(what = RSLD) |> expect_error()
  c(ARGS[-3L], TRUE)  |> do.call(what = RSLD) |> expect_no_error()
  c(ERR[-3L]) |> do.call(what = RIDGESearchLambdaDichomotomy) |> expect_error()
  RIDGESearchLambdaDichomotomy(ARGS[[1]], ARGS[[2]]) |> expect_no_error()
  RIDGESearchLambdaDichomotomy(ERR[[1]], ERR[[2]]) |> expect_error()
})







