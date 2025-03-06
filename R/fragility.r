#' Compute the normalized fragility row for adjacency matrix A
#'
#' @param A Numeric. Adjacency Matrix
#' @param nSearch Integer. Number of eigenvalues tried to find the minimum norm vector
#' @param normalize Logical. If TRUE, the fragility row is normalized
fragilityRow <- function(A, nSearch = 100, normalize = TRUE) {
  nel <- ncol(A)
  ek <- diag(nel)
  b <- matrix(c(0, -1), 2L)
  fragNorm <- rep(0L, nel)
  omvec <- seq(0L, 1L, length.out = nSearch + 1L)[-1L]
  lambdas  <- sqrt(1 - omvec^2) + omvec * 1i
  ## (A - (sigma + j * omega)*I)^-1
  iMats <- lapply(lambdas, \(L) t(solve(A - L * ek)))
  for (i in seq_len(nel)) {
    minNorm <- 100000
    item <- ek[i, , drop = FALSE]
    for (k in seq_len(nSearch)) {
      argument <- item %*% iMats[[k]]
      B <- rbind(Im(argument), Re(argument))
      ## B^T * (B * B^T)^-1 * b
      prov <- t(B) %*% solve(B %*% t(B)) %*% b
      provNorm <- norm(prov, type = "2")
      if (provNorm < minNorm) minNorm <- provNorm
    }
    fragNorm[i] <- minNorm
  }
  if (!normalize) return(fragNorm)
  maxf <- max(fragNorm)
  (maxf - fragNorm) / maxf
}

#' Compute quantiles, mean and standard deviation for two electrodes group marked as soz non marked as soz
#'
#' @param frag Matrix or Fragility object. Either a matrix with row as Electrode names and Column as fragility index, or a Fragility object from \code{calcAdjFrag}

#' @param sozID Integer.  Vector soz electrodes (for good electrodes)
#'
#'
#' @return list of 5 items with quantile matrix, mean and sdv from both electrodes groups
#' @export
#'
#' @examples
#' data("pt01Fragm1sp2s")
#' data("pt01Epochm1sp2s")
#' sozindex<-attr(pt01Epochm1sp2s,"sozindex")
#' pt01fragstat<-fragStat(frag=pt01Fragm1sp2s, sozID=sozindex)
fragStat <- function(frag, sozID) {
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


predictRidge <- function(xt, A) {
    ## the data matrix
    if (nrow(A) == ncol(A) + 1) {
        x <- cbind(1, as.matrix(xt))
    } else {
        x <- as.matrix(xt)
    }
    x %*% A
}
