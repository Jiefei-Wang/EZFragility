# I rewrote ridge script for clarity and efficiency
# For clarity, it is at least more clear to me.
# For efficiency, it definitely runs faster.
RIDGE <- function(x, x1, l) {
  n <- ncol(x)
  cb <- \(i) { 
    glmnet::glmnet(
      x = x, y = x1[, i], lambda = l, nlambda = 1, alpha = 0, 
      standardize = FALSE, intercept = FALSE, family = "gaussian", 
      type.gaussian = "covariance"
    )[["beta"]]@x
  }
  A <- vapply(seq_len(n), cb, numeric(n), USE.NAMES = FALSE)
  isStable <- (eigen(A, FALSE, TRUE)$values |> Mod() |> max()) < 1.0
  structure(A, lambda = l, stable = isStable)
}

RIDGESearchLambdaDichomotomy <- function(xt, xtp1) {
  if (!identical(dim(xt), dim(xtp1))) stop("Unmatched dimension")
  low = lmbd <- 1e-4
  high <- low * 2^3
  A <- RIDGE(xt, xtp1, lmbd);
  if (!attr(A, "stable")) {
    A <- RIDGE(xt, xtp1, high)
    while (!attr(A, "stable")) {
      low <- high
      high <- high * 2^2
      A <- RIDGE(xt, xtp1, high)
    }
    for (i in seq_len(20L)) {
      lmbd <- (low + high) * .5
      A <- RIDGE(xt, xtp1, lmbd)
      if (attr(A, "stable")) {
        high <- lmbd
        if (.5 * (1 - low/high) < .2) break
      }
      else low <- lmbd
    }
  }
  structure(A, lambda = lmbd)
}

RidgeRSQ <- function(xt, xtp1, A) {
  E <- xtp1 - xt %*% A
  yc <- t(t(xtp1) - colMeans(xtp1))
  sse <- diag(crossprod(E))
  sst <- diag(crossprod(yc))
  1.0 - sse / sst
}
