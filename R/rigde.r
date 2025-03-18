#' fit a generalized linear model to compute adjacency matrix A
#'
#' A x(t) = x(t+1)
#'
#' @param xt matrix. iEEG time series for a given window,
#' with time points as rows and electrodes names as columns
#' @param xtp1 matrix. the iEEG time serie at the next time point,
#' with time points as rows and electrodes names as columns
#' @param lambda Numeric Vector. A user supplied lambda sequence.
#' @return adjacency matrix A
ridge <- function(xt, xtp1, lambda) {
  Xsv <- svd(xt)
  d <- Xsv$d
  uT <- t(Xsv$u)
  v <- Xsv$v

  n <- nrow(xt)
  p <- ncol(xt)
  lmbd <- n * lambda

  .cback <- function(i) {
    y <- xtp1[, i]
    Lscaled <- lmbd * sum(y^2 / n)^-.5
    dw <- d * (d^2 + Lscaled)^-1
    drop(v %*% (dw * uT %*% y))
  }

  A <- vapply(seq_len(p), .cback, numeric(p), USE.NAMES = FALSE)
  isStable <- (eigen(A, FALSE, TRUE)$values |> Mod() |> max()) < 1.0
  structure(A, lambda = lambda, stable = isStable)
}


#' computes R2
#'
#' @inheritParams ridge
#' @param A adjacency matrix
#'
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
#'
#' @inheritParams ridge
#'
#' @return adjacency matrix Afin with lambda as attribute
ridgeSearch <- function(xt, xtp1, lambda = NULL) {
  if (!identical(dim(xt), dim(xtp1))) stop("Unmatched dimension")
  if (!is.null(lambda)) {
    A <- ridge(xt, xtp1, lambda)
    return(structure(A, lambda = lambda))
  }
  low <- lambda <- 1e-4
  high <- 10
  A <- ridge(xt, xtp1, lambda)
  if (!attr(A, "stable")) {
    for (i in seq_len(20L)) {
      l <- (low + high) * .5
      A_tmp <- ridge(xt, xtp1, l)
      if (attr(A_tmp, "stable")) {
        high <- lambda <- l
        A <- A_tmp
      } else {
        low <- l
      }
    }
  }
  structure(A, lambda = lambda)
}
