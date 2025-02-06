# Same are ridge script, I rewrote it for clarity and efficiency
# There is not need for two fragilityRow functions.
# It can be one function with normalize as boolean argument
FragilityROW <- function(A, nSearch = 100, normalize = TRUE) {
  nel <- ncol(A)
  ek <- diag(nel)
  b <- matrix(c(0, -1), 2L)
  fragNorm <- rep(0L, nel)
  omvec <- seq(0L, 1L, length.out = nSearch + 1L)[-1L]
  lambdas  <- sqrt(1 - omvec^2) + omvec * 1i
  iMats <- lapply(lambdas, \(L) t(solve(A - L * ek)))
  for (i in seq_len(nel)) {
    minNorm <- 100000
    item <- ek[i, , drop = FALSE]
    for (k in seq_len(nSearch)) {
      argument <- item %*% iMats[[k]]
      B <- rbind(Im(argument), Re(argument))
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

# I rewrote this one just for clarity
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
