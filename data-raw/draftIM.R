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


## We will use data from -1s to 2s around the seizure onset

#
# si   <- i0 + (iw - 1L) * t_step
#
# fn <- \(win, tst, dat, mapDatRangeInSeconds) {
#   self <- environment()
#   dat <- dat;
#   w <- win;
#   u <- tst;
#   n <- nrow(dat);
#   tr <- mapDatRangeInSeconds
#   tm <- do.call(seq, as.list(c(tr, by = 1e-3)))[-(n + 1L)]
#   sb <- \(x, y = NULL) {
#     if (length(x) == 2L) {
#       y <- x[2L]
#       x <- x[1L]
#     }
#     else if (is.null(y)) return(dat[x  <= tm, ])
#     else if (is.null(x)) return(dat[tm <  y,  ])
#     dat[x <= tm & tm < y, ]
#   }
#   Step2Idx <- \(i) w + u * (i - 1L)
#   self
# }
#
# SubsetFragilityObj <- \(win, tst, frag, mapDatRangeInSeconds) {
#   self <- environment()
#   frag <- frag;
#   w <- win;
#   u <- tst;
#   dat <- frag@ieegts
#   n <- nrow(dat);
#   tr <- mapDatRangeInSeconds
#   steps <- dim(frag@adj)[3L]
#   Step2Idx <- \(i) w + u * (i - 1L)
#   idx2Step <- \(idx) (idx - w) %/% u + 1L
#   tm <- do.call(seq, as.list(c(tr, by = 1e-3)))[-(n + 1L)]
#   sb <- \(x, y = NULL) {
#     if (length(x) == 2L) {
#       y <- x[2L]
#       x <- x[1L]
#     }
#     else if (is.null(y)) return(dat[x  <= tm, ])
#     else if (is.null(x)) return(dat[tm <  y,  ])
#     dat[x <= tm & tm < y, ]
#   }
#   self
# }
#
# ff <- \(win, tst, frag, mapDatRangeInSeconds) {
#   self <- environment()
#   frag <- frag;
#   w <- win;
#   u <- tst;
#   n <- nrow(dat);
#   tr <- mapDatRangeInSeconds
#   tm <- do.call(seq, as.list(c(tr, by = 1e-3)))[-(n + 1L)]
#   sb <- \(x, y = NULL) {
#     if (length(x) == 2L) {
#       y <- x[2L]
#       x <- x[1L]
#     }
#     else if (is.null(y)) return(dat[x  <= tm, ])
#     else if (is.null(x)) return(dat[tm <  y,  ])
#     dat[x <= tm & tm < y, ]
#   }
#   Step2Idx <- \(i) w + u * (i - 1L)
#   self
# }
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# FragSlotOnly <- pt01Frag
#
# exList <- slotNames(FragSlotOnly) |>
#   (\(x) x[x != "frag"])() |>
#   lapply(\(i) substitute(FragSlotOnly@x <- NULL, list(x = i)))
#
# for (i in seq_along(exList)) eval(exList[[i]])
#
#
#
# lapply(sl, \(i) substitute(frag1@x <- NULL, list(x = i)))
#
# argL <- list(
#   ieegts = pt01Epoch[9001:12000,],
#   t_window = 250,
#   t_step = 125,
#   lambda = NULL,
#   nSearch = 100
# )
#
# system.time({ fragNew <- do.call(GetAdjFrag, argL) })
# system.time({ fragOld <- do.call(calc_adj_frag, argL) })
#
# rbenchmark::benchmark(
#   calc_adj_frag = do.call(calc_adj_frag, argL),
#   GetAdjFrag = do.call(GetAdjFrag, argL),
#   replications = 1
# )
#
# res = do.call(GetAdjFrag, argL)
#
# scaling <- 10^floor(log10(max(iEEGts)))
# X  <- iEEGts / scaling
#
# objf <- function(A, lambda, xt, xtp1) {
#   n <- nrow(xt)
#   return( sum((xtp1 - xt %*% A)^2) + lambda * sum(A^2) )
# }
#
# X <- pt01Epoch[f(5), ]
# Y <- pt01Epoch[f(5) + 1, ]
#
#
#
# R1 <- function(xt, xtp1) {
#   if (!identical(dim(xt), dim(xtp1))) stop("Unmatched dimension")
#   e <- environment()
#   low = current = lambdaopt <- 0.0001
#   high = 10
#   A <- RIDGE(xt, xtp1, current)
#   isStable <- \(e) ( eigen(e$A)$values |> Mod() |> max() ) < 1
#   updateRange <- \(e) {
#     if (isStable(e)) {
#       e$high <- e$current
#       e$lambdaopt <- e$current
#     }
#     else {
#       e$low <- e$current
#       e$lambdaopt <- e$high
#     }
#     (e$high - e$low) < .01
#   }
#   if (!isStable(e)) {
#     for (i in seq_len(20L)) {
#       current <- (low + high) * .5
#       A <- RIDGE(xt, xtp1, current)
#       if (updateRange(e)) break
#     }
#   }
#   if (current != lambdaopt) A <- RIDGE(xt, xtp1, lambdaopt)
#   structure(A, lambda = lambdaopt)
# }
#
#
#
# R2 <- function(xt, xtp1) {
#   if (!identical(dim(xt), dim(xtp1))) stop("Unmatched dimension")
#   e <- environment()
#   low = current = lambdaopt <- 0.0001
#   high = 10
#   A <- RIDGE(xt, xtp1, current)
#
#   isStable <- \() ( eigen(e$A)$values |> Mod() |> max() ) < 1
#   updateRange <- \() {
#     if (isStable()) {
#       e$high <- e$current
#       e$lambdaopt <- e$current
#     }
#     else {
#       e$low <- e$current
#       e$lambdaopt <- e$high
#     }
#     (e$high - e$low) < .01
#   }
#
#   if (!isStable()) {
#     for (i in seq_len(20L)) {
#       current <- (low + high) * .5
#       A <- RIDGE(xt, xtp1, current)
#       if (updateRange()) break
#     }
#   }
#   if (current != lambdaopt) A <- RIDGE(xt, xtp1, lambdaopt)
#   structure(A, lambda = lambdaopt)
# }
#
# ieegts <- pt01Epoch[1:750, ]
# # ieegts <- ieegts/ 10^floor(log10(max(ieegts)))
# argL <- list(
#   ieegts = ieegts,
#   t_window = 250,
#   t_step = 125,
#   lambda = NULL,
#   nSearch = 1
# )
#
# fragNew <- do.call(GetAdjFrag, argL)
#
# f <- \(i) seq_len(249) + -125 + 125 * i
#
# isst <- \(A) ( eigen(A)$values |> Mod() |> max() ) < 1
#
#
# ff <- \(i) {
#   x <- ieegts[f(i), ]
#   y <- ieegts[f(i) + 1, ]
#   A <- RIDGE(x, y, 0.5)
#   eigen(A, FALSE, TRUE)$values |> Mod() |> max()
# }
#
# sapply(1:5, ff)
#
# ff <- \(i) {
#   x <- ieegts[f(i), ]
#   y <- ieegts[f(i) + 1, ]
#   attr(RIDGESearchLambdaDichomotomy(x, y), "lambda")
# }
#
#
#
# i <- 4
# x <- ieegts[f(i), ]
# y <- ieegts[f(i) + 1, ]
# L <- 0.5
#
#
# A <- RIDGE(x, y, L)
# (Fx0 <- eigen(A)$values |> Mod() |> max())
#
# dL <- 1e-8
# B <- RIDGE(x, y, L + dL)
# (Fx1 <- eigen(B)$values |> Mod() |> max())
# (dFx <- (Fx1 - Fx0) / dL)
# (L <- L - 0.1 * Fx0/dFx)
#
# x <- ieegts[f(1), ]
# y <- ieegts[f(1) + 1, ]
#
# ll <- diag(84) * 1e-4
#
# solve(t(x) %*% x + ll) %*% t(x) %*% y[,1]

#
# ridge1 <- function(x, y, l) {
#   n <- ncol(x)
#   nc <- ncol(y)
#   rhs <- solve(t(x) %*% x + diag(n) * l) %*% t(x)
#   cb <- \(i) drop(rhs %*% y[, i])
#   vapply(seq_len(nc), cb, numeric(n), USE.NAMES = FALSE) |> structure(lambda = l)
# }
#
# myRidge <- function(x, y, l) {
#   n <- ncol(x)
#   cb <- \(i) c(solve(t(x) %*% x + diag(n) * l) %*% t(x) %*% y[, i])
#   vapply(seq_len(n), cb, numeric(n), USE.NAMES = FALSE) |> structure(lambda = l)
# }
#
# y <- xtp1[, 1, drop = FALSE]
# y_mean <- mean(y)
# y_sd <- sqrt(sum((y - y_mean)^2) / length(y))
# # y_sd <- sqrt(sum(y^2) /length(y))
# y <- y  / y_sd

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














