## We will use data from -1s to 2s around the seizure onset


argL <- list(
  ieegts = pt01Epoch[9001:12000,],
  t_window = 250,
  t_step = 125,
  lambda = NULL,
  nSearch = 100
)

system.time({ fragNew <- do.call(GetAdjFrag, argL) })
system.time({ fragOld <- do.call(calc_adj_frag, argL) })

rbenchmark::benchmark(
  calc_adj_frag = do.call(calc_adj_frag, argL),
  GetAdjFrag = do.call(GetAdjFrag, argL),
  replications = 1
)

res = do.call(GetAdjFrag, argL)

scaling <- 10^floor(log10(max(iEEGts)))
X  <- iEEGts / scaling

objf <- function(A, lambda, xt, xtp1) {
  n <- nrow(xt)
  return( sum((xtp1 - xt %*% A)^2) + lambda * sum(A^2) )
}

X <- pt01Epoch[f(5), ]
Y <- pt01Epoch[f(5) + 1, ]



R1 <- function(xt, xtp1) {
  if (!identical(dim(xt), dim(xtp1))) stop("Unmatched dimension")
  e <- environment()
  low = current = lambdaopt <- 0.0001
  high = 10
  A <- RIDGE(xt, xtp1, current)
  isStable <- \(e) ( eigen(e$A)$values |> Mod() |> max() ) < 1
  updateRange <- \(e) {
    if (isStable(e)) {
      e$high <- e$current
      e$lambdaopt <- e$current
    }
    else {
      e$low <- e$current
      e$lambdaopt <- e$high
    }
    (e$high - e$low) < .01
  }
  if (!isStable(e)) {
    for (i in seq_len(20L)) {
      current <- (low + high) * .5
      A <- RIDGE(xt, xtp1, current)
      if (updateRange(e)) break
    }
  }
  if (current != lambdaopt) A <- RIDGE(xt, xtp1, lambdaopt)
  structure(A, lambda = lambdaopt)
}



R2 <- function(xt, xtp1) {
  if (!identical(dim(xt), dim(xtp1))) stop("Unmatched dimension")
  e <- environment()
  low = current = lambdaopt <- 0.0001
  high = 10
  A <- RIDGE(xt, xtp1, current)

  isStable <- \() ( eigen(e$A)$values |> Mod() |> max() ) < 1
  updateRange <- \() {
    if (isStable()) {
      e$high <- e$current
      e$lambdaopt <- e$current
    }
    else {
      e$low <- e$current
      e$lambdaopt <- e$high
    }
    (e$high - e$low) < .01
  }

  if (!isStable()) {
    for (i in seq_len(20L)) {
      current <- (low + high) * .5
      A <- RIDGE(xt, xtp1, current)
      if (updateRange()) break
    }
  }
  if (current != lambdaopt) A <- RIDGE(xt, xtp1, lambdaopt)
  structure(A, lambda = lambdaopt)
}

ieegts <- pt01Epoch[1:750, ]
# ieegts <- ieegts/ 10^floor(log10(max(ieegts)))
argL <- list(
  ieegts = ieegts,
  t_window = 250,
  t_step = 125,
  lambda = NULL,
  nSearch = 1
)

fragNew <- do.call(GetAdjFrag, argL)

f <- \(i) seq_len(249) + -125 + 125 * i

isst <- \(A) ( eigen(A)$values |> Mod() |> max() ) < 1


ff <- \(i) {
  x <- ieegts[f(i), ]
  y <- ieegts[f(i) + 1, ]
  A <- RIDGE(x, y, 0.5)
  eigen(A, FALSE, TRUE)$values |> Mod() |> max()
}

sapply(1:5, ff)

ff <- \(i) {
  x <- ieegts[f(i), ]
  y <- ieegts[f(i) + 1, ]
  attr(RIDGESearchLambdaDichomotomy(x, y), "lambda")
}



i <- 4
x <- ieegts[f(i), ]
y <- ieegts[f(i) + 1, ]
L <- 0.5


A <- RIDGE(x, y, L)
(Fx0 <- eigen(A)$values |> Mod() |> max())

dL <- 1e-8
B <- RIDGE(x, y, L + dL)
(Fx1 <- eigen(B)$values |> Mod() |> max())
(dFx <- (Fx1 - Fx0) / dL)
(L <- L - 0.1 * Fx0/dFx)

x <- ieegts[f(1), ]
y <- ieegts[f(1) + 1, ]

ll <- diag(84) * 1e-4

solve(t(x) %*% x + ll) %*% t(x) %*% y[,1]



myRigde <- function(x, x1, l) {
  n <- ncol(x)
  cb <- \(i) c(solve(t(x) %*% x + diag(n) * l) %*% t(x) %*% y[,i])
  vapply(seq_len(n), cb, numeric(n), USE.NAMES = FALSE) |> structure(lambda = l)
}

#' Compute quantiles, mean and standard deviation for two electrodes group marked as soz non marked as soz
#'
#' @param frag Matrix or Fragility object. Either a matrix with row as Electrode names and Column as fragility index, or a Fragility object from \code{calc_adj_frag}
#' @param elecsoz Integer.  Vector soz electrodes (for good electrodes)
#'
#'
#' @return list of 5 items with quantile matrix, mean and sdv from both electrodes groups
#' @export
#'
#' @examples
#' data("pt01Fragility")
#' data("elecsoz")
#' fragstat<-frag_stat(frag=pt01Fragility, elecsoz=elecsoz)

frag_stat <- function(frag, sozID) {
  if (is(frag, "Fragility")) frag <- frag$frag
  if (!inherits(frag, "matrix")) stop("Frag must be matrix or Fragility object")
  steps <- ncol(frag)
  sozCID <- which(!(seq_len(nrow(frag)) %in% sozID)) 
  hmapSOZ  <- frag[sozID,  , drop = FALSE]
  hmapSOZC <- frag[sozCID, , drop = FALSE]
  muSOZ  <- colMeans(hmapSOZ)
  muSOZC <- colMeans(hmapSOZC)
  sdSOZ  <- apply(hmapSOZ,  2L, sd)
  sdSOZC <- apply(hmapSOZC, 2L, sd)
  Q <- seq(.1, 1, by = .1)
  qmatrix <- rbind(
    apply(hmapSOZ,  2, quantile, Q),
    apply(hmapSOZC, 2, quantile, Q)
  )
  rowPrefix <- rep(c("SOZ", "SOZC"), each = 10)
  dimN <- dimnames(qmatrix)
  dimnames(qmatrix) <- list(
    Quantiles = paste0(rowPrefix, dimN[[1L]]),
    Step      = dimN[[2L]]
  )
  FragStat(
    qmatrix   = qmatrix,
    cmeansoz  = muSOZ,
    cmeansozc = muSOZC,
    csdsoz    = sdSOZ,
    csdsozc   = sdSOZC
  )
}
