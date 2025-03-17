fragIdObj <- \(ieegts, mframe, period, origin = 0) {
  self <- environment()
  n = nrow(ieegts)
  endTime <- n * 1e-3 + origin
  mapStep2Row <- \(s) mframe + period * (s - 1L)
  mapRow2Step <- \(i) (i - 1L) %/% period + 1 * (i < mframe / 2)
  init <- \(t0 = origin) {
    self$origin <- t0
    self$endTime <- n * 1e-3 + origin
    self$mapTime2Row <- \(x) (round(x, 3) - origin) * 1e3
    self$getIds <- \(x = origin, y = endTime) {
      if (length(x) == 2L) {y <- x[2L]; x <- x[1L]}
      stopifnot(x < y, origin <= x, y <= endTime)
      a <- mapTime2Row(x) + 1
      z <- mapTime2Row(y)
      list(rows = a:z, steps = mapRow2Step(a):mapRow2Step(z))
    }
    invisible(self)
  }
  setOrigin <- \(newOrigin) init(newOrigin)
  init()
}

def <- \(x) {
  l <- x
  \(...) {
    inp <- list(...)
    unNamed <- list()
    for (i in seq_along(inp)) {
      n <- names(inp[i])
      x <- if (is.symbol(inp[[i]])) eval(inp[[i]]) else inp[[i]]
      if (is.atomic(x)) x <- list(x)
      if (length(nchar(n))) l[n] <- x
      else unNamed <- c(unNamed, if (is.atomic(x)) list(x) else x)
    }
    c(unNamed, l)
  }
}

glmNetLambda <- \(y, inputLambda) {
  n <- length(y)
  n * inputLambda * sum(y^2 / n)^-.5
}

rdgUnit <- \(x, y, l) {
  lmat <- diag(ncol(x)) * l
  beta <- solve(t(x) %*% x + lmat) %*% t(x) %*% y
  drop(beta)
}

ridge <- function(xt, xtp1, lambda) {
  if (!inherits(xtp1, "matrix")) xtp1 <- as.matrix(xtp1)
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

  A <- vapply(seq_len(ncol(xtp1)), .cback, numeric(p), USE.NAMES = FALSE)
  isStable <- logical()
  isStable <-  if (nrow(A) == ncol(A)) (eigen(A, FALSE, TRUE)$values |> Mod() |> max()) < 1.0
  # isStable <- (eigen(A, FALSE, TRUE)$values |> Mod() |> max()) < 1.0
  structure(A, lambda = lambda, stable = isStable)
}

ridge1 <- function(xt, xtp1, lambda) {
  ei <- eigen(t(xt) %*% xt)
  u <- ei$vectors
  d <- ei$values
  rhs <- t(u) %*% t(xt)


  n <- nrow(xt)
  p <- ncol(xt)
  lmbd <- n * lambda

  .cback <- function(i) {
    y <- xtp1[, i]
    Lscaled <- lmbd * sum(y^2 / n)^-.5
    dw <- (d + Lscaled)^-1
    drop(u %*% (dw * rhs %*% y))
  }

  A <- vapply(seq_len(p), .cback, numeric(p), USE.NAMES = FALSE)
  isStable <- (eigen(A, FALSE, TRUE)$values |> Mod() |> max()) < 1.0
  structure(A, lambda = lambda, stable = isStable)
}


require(glmnet)
data("pt01Epochm1sp2s")
X <-  pt01Epochm1sp2s / (10^floor(log10(max(pt01Epochm1sp2s))))
# X <-  pt01Epochm1sp2s / sd(pt01Epochm1sp2s)
i = 0
lmbd <- 1e-4

x = xt <- X[1:249 + 125 * i, ]
xtp1   <- X[2:250 + 125 * i, ]
y      <- xtp1[, 1]

dargs <- list(
  x = xt,
  y = xtp1[, 1],
  lambda = lmbd,
  alpha = 0,
  thresh = 1e-15,
  maxit = 1e6,
  standardize = FALSE,
  intercept = FALSE
)

vL <- def(dargs)

{
  fit <- do.call(glmnet, vL(thresh = 1e-25));
  myBeta <- ridge(xt, y, lmbd)
  cbind(fit$beta@x, myBeta)
  all.equal(fit$beta@x, vapply(myBeta, I, .0))
}














