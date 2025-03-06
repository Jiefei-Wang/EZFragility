# setup ------------------------------------------------------------------------
set.seed(123L)
mat <- matrix(rnorm(1e3L), ncol = 10L)
ARGS = ERR <- list(xt = mat[1L:4L, ], xtp1 = mat[1L:4L + 1L, ], 0.1)
ERR$xtp1 <- ERR$xtp1[, -1L]

# fragilityRow -------------------------------------------------------
test_that("fragilityRow", {
  input <- do.call(ridge, ARGS)
  fragilityRow(input) |> expect_no_error()
})

# ridge ------------------------------------------------------------------------
test_that("ridge/ridgeR2", {
  do.call(ridge, ERR)  |> expect_error()
  A <- do.call(ridge, ARGS) |> expect_no_error()
  do.call(ridgeR2, c(ARGS[-3L], list(A))) |> expect_no_error()
})

# ridgesearchlambdadichomotomy -------------------------------------------------
test_that("ridgesearchlambdadichomotomy", {
  do.call(ridgesearchlambdadichomotomy, c(ERR[-3L], FALSE)) |> expect_error()
  do.call(ridgesearchlambdadichomotomy, c(ARGS[-3L], TRUE)) |> expect_no_error()
})


# internal data -------------------------------------------------
test_that("Data consistency across versions", {
  set.seed(1)
  data <- matrix(rnorm(800), nrow = 40)
  t_window <- 10
  t_step <- 5
  lambda <- 0.1
  frag <- calcAdjFrag(ieegts = data, window = t_window, step = t_step, lambda = lambda)
  expect_equal(frag, testFrag)
})
