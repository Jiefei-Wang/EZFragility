#' Compute the normalized fragility row for adjacency matrix A
#' 
#' The matrix A is used for the regression: A * x(t) = x(t+1)
#'
#' @param A Numeric. Adjacency Matrix
#' @param nSearch Integer. Number of eigenvalues tried to find the minimum norm vector
#' @param normalize Logical. If TRUE, the fragility row is normalized
fragilityRow <- function(A, nSearch = 100, normalize = TRUE) {
    elecNum <- nrow(A)
    ek <- diag(elecNum)
    b <- matrix(c(0, -1), 2L)
    fragNorm <- rep(0L, elecNum)
    omvec <- seq(0L, 1L, length.out = nSearch + 1L)[-1L]
    lambdas <- sqrt(1 - omvec^2) + omvec * 1i
    ## (A - (sigma + j * omega)*I)^-1
    iMats <- lapply(lambdas, \(L) solve(A - L * ek))
    for (i in seq_len(elecNum)) {
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
    if (!normalize) {
        return(fragNorm)
    }
    maxf <- max(fragNorm)
    ## TODO: find another way to normalize the fragility row
    ## The current implementation reverse the meaning of the fragility
    ## TODO: add both frag and fragNormalized to Fragility
    (maxf - fragNorm) / maxf
}

#' Compute quantiles, mean and standard deviation for two electrodes group 
#'
#' @param frag Matrix or Fragility object. Either a matrix with row as Electrode names and Column as fragility index, or a Fragility object from \code{calcAdjFrag}

#' @param groupIndex Integer.  Vector group electrodes (for good electrodes)
#'
#'
#' @return list of 5 items with quantile matrix, mean and sdv from both electrodes groups
#'
#' @examples
#' data("pt01Frag")
#' data("pt01EcoG")    
#' ## sozNames is the name of the electrodes we assume are in the SOZ
#' sozNames <- metaData(pt01EcoG)$sozNames
#' pt01fragstat <- fragStat(frag = pt01Frag, groupIndex = sozNames)
#' @export 
fragStat <- function(frag, groupIndex,groupName="SOZ") {
## TODO: support grouped and ungrouped fragility statistics (Not now, but for the future)
    if (is(frag, "Fragility")) frag <- frag$frag
    if (!inherits(frag, "matrix")) stop("Frag must be matrix or Fragility object")
    steps <- ncol(frag)
    groupCID <- which(!(seq_len(nrow(frag)) %in% groupIndex))
    hmapGroup <- frag[groupIndex, , drop = FALSE]
    hmapREF <- frag[groupCID, , drop = FALSE]
    meanGroup <- colMeans(hmapGroup)
    meanRef <- colMeans(hmapREF)
    sdGroup <- apply(hmapGroup, 2L, sd)
    sdRef <- apply(hmapREF, 2L, sd)
    Q <- seq(.1, 1, by = .1)
    qmatrix <- rbind(
        apply(hmapGroup, 2, quantile, Q),
        apply(hmapREF, 2, quantile, Q)
    )
    
    rowPrefix <- rep(c(groupName, "REF"), each = 10)
    dimN <- dimnames(qmatrix)
    dimnames(qmatrix) <- list(
        Quantiles = paste0(rowPrefix, dimN[[1L]]),
        Step      = dimN[[2L]]
    )
    FragStat(
        qmatrix   = qmatrix,
        meanGroup  = meanGroup,
        meanRef = meanRef,
        sdGroup    = sdGroup,
        sdRef   = sdRef
    )
}

