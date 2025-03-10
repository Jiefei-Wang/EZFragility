#' fit a generalized linear model to compute adjacency matrix A
#'
#' A x(t) = x(t+1)
#' @param xt matrix. iEEG time series for a given window,
#' with time points as rows and electrodes names as columns
#' @param xtp1 matrix. the iEEG time serie at the next time point,
#' with time points as rows and electrodes names as columns
#' @param lambda Numeric Vector. A user supplied lambda sequence.
#' @return adjacency matrix A
ridge <- function(xt, xtp1, lambda) {
  n <- ncol(xt)
  cb <- \(i) {
    glmnet::glmnet(x = xt, y = xtp1[, i], lambda = lambda, alpha = 0,
      standardize = FALSE, intercept = FALSE)[["beta"]]@x
  }
  A <- vapply(seq_len(n), cb, numeric(n), USE.NAMES = FALSE)
  isStable <- (eigen(A, FALSE, TRUE)$values |> Mod() |> max()) < 1.0
  structure(A, lambda = lambda, stable = isStable)
}
#' computes R2
#'
#' @inheritParams ridge
#' @param A adjacency matrix
ridgeR2 <- function(xt, xtp1, A) {
  nel <- ncol(xt)
  ypredMat <- predictRidge(xt, A)

  R2 <- rep(0, nel)
  for (i in seq_len(nel)) {
    y <- xtp1[, i]
    ypred <- ypredMat[, i]
    sst <- sum((y - mean(y))^2)
    sse <- sum((ypred - y)^2)
    rsq <- 1 - sse / sst
    R2[i] <- rsq
  }
  R2
}
#' Ridge Regression for Electrode Readings
#'
#' Ridge regression to compute matrix adjancency matrix A such as A xt = xtpt1
#' the lambda parmeter is found by dichotomy such that A is stable
#' (all eigenvalues have a norm less than one)
#' @inheritParams ridge
#' @return adjacency matrix Afin with lambda as attribute
RidgeSearchLambdaDichomotomy <- function(xt, xtp1) {
  if (!identical(dim(xt), dim(xtp1))) stop("Unmatched dimension")
  low = lambda <- 1e-4
  high <- 10
  A <- ridge1(xt, xtp1, lambda);
  if (!attr(A, "stable")) {
    for (i in seq_len(20L)) {
      l <- (low + high) * .5
      A_tmp <- ridge(xt, xtp1, l)
      if (attr(A_tmp, "stable")) {
        high = lambda <- l
        A <- A_tmp
      }
      else low <- l
    }
  }
  structure(A, lambda = lambda)
}
