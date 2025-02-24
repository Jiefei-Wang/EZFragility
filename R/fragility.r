#' Compute the normalized fragility row for adjacency matrix A
#' 
#' @param A Numeric. Adjacency Matrix  
#' @param n_search Integer. Number of eigenvalues tried to find the minimum norm vector 
fragilityRowNormalized <- function(A, n_search = 100) {
  ## The adjacency matrix A here is a transpose of the
  ## adjacency matrix in the original paper
  nel <- ncol(A)
  e <- Mod(eigen(A)$values)
  me <- max(e)

  if (me >= 1) {
    #return(0)
  }


  fragcol <- matrix(0, nel, nel)
  fragNorm <- rep(0, nel)
  omvec <- seq(0, 1, length.out = n_search + 1)[-1]

  b <- c(0, -1)
  ## for each electrode
  for (i in 1:nel) {
    ## indicate which electrode is disturbed (ith)
    ek <- rep(0, nel)
    ek[i] <- 1
    tek <- t(ek)
    minNorm <- 100000
    minPerturbColumn <- NA
    for (k in seq_len(n_search)) {
      ## imaginary part
      om <- omvec[k]
      ## real part
      sigma <- sqrt(1 - om^2)
      ## target eigenvalue
      lambda <- complex(real = sigma, imaginary = om)
      ## A - (sigma + j* omega)*I
      mat <- A - lambda * diag(nel)
      imat <- t(solve(mat))

      argument <- tek %*% imat
      B <- rbind(Im(argument), Re(argument))
      ## B^T*(B*B^T)^-1*b
      invBtB <- solve(B %*% t(B))
      prov <- t(B) %*% invBtB %*% b

      sigma_hat <- ek %*% t(prov)

      ## validation
      if (FALSE) {
        A2 <- A + sigma_hat
        e2 <- eigen(A2)$values
        closestIndex <- which.min(abs(e2 - lambda))
        e2[closestIndex]
      }

      norm_sigma_hat <- norm(sigma_hat, type = "2")
      if (norm_sigma_hat < minNorm) {
        minPerturbColumn <- prov
        minNorm <- norm(prov, type = "2")
      }
    }

    fragcol[, i] <- minPerturbColumn
    fragNorm[i] <- minNorm
  }

  maxf <- max(fragNorm)
  fragNorm2 <- (maxf - fragNorm) / maxf

  return(fragNorm2)
}

#' Compute the fragility row for adjacency matrix A
#'
#' @inheritParams fragilityRowNormalized
fragilityRow <- function(A, n_search = 100) {
  ## The adjacency matrix A here is a transpose of the
  ## adjacency matrix in the original paper
  nel <- ncol(A)
  e <- Mod(eigen(A)$values)
  me <- max(e)

  if (me >= 1) {
    return(0)
  }


  fragcol <- matrix(0, nel, nel)
  fragNorm <- rep(0, nel)
  omvec <- seq(0, 1, length.out = n_search + 1)[-1]

  b <- c(0, -1)
  ## for each electrode
  for (i in 1:nel) {
    ## indicate which electrode is disturbed (ith)
    ek <- rep(0, nel)
    ek[i] <- 1
    tek <- t(ek)
    minNorm <- 100000
    minPerturbColumn <- NA
    for (k in seq_len(n_search)) {
      ## imaginary part
      om <- omvec[k]
      ## real part
      sigma <- sqrt(1 - om^2)
      ## target eigenvalue
      lambda <- complex(real = sigma, imaginary = om)
      ## A - (sigma + j* omega)*I
      mat <- A - lambda * diag(nel)
      imat <- t(solve(mat))

      argument <- tek %*% imat
      B <- rbind(Im(argument), Re(argument))
      ## B^T*(B*B^T)^-1*b
      invBtB <- solve(B %*% t(B))
      prov <- t(B) %*% invBtB %*% b

      sigma_hat <- ek %*% t(prov)

      ## validation
      if (FALSE) {
        A2 <- A + sigma_hat
        e2 <- eigen(A2)$values
        closestIndex <- which.min(abs(e2 - lambda))
        e2[closestIndex]
      }

      norm_sigma_hat <- norm(sigma_hat, type = "2")
      if (norm_sigma_hat < minNorm) {
        minPerturbColumn <- prov
        minNorm <- norm(prov, type = "2")
      }
    }

    fragcol[, i] <- minPerturbColumn
    fragNorm[i] <- minNorm
  }

  return(fragNorm)
}

#' Compute quantiles, mean and standard deviation for two electrodes group marked as soz non marked as soz
#'
#' @param frag Matrix or Fragility object. Either a matrix with row as Electrode names and Column as fragility index, or a Fragility object from \code{calc_adj_frag}
#' @param sozID Integer.  Vector soz electrodes (for good electrodes)
#' 
#'
#' @return list of 5 items with quantile matrix, mean and sdv from both electrodes groups
#' @export
#'
#' @examples
#' data("pt01Frag")
#' data("pt01Epoch")
#' sozindex<-attr(pt01Epoch,"sozindex")
#' fragstat<-frag_stat(frag=pt01Frag, sozID=sozindex)
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

# Compute quantiles, mean and standard deviation for two electrodes group marked as soz non marked as soz
# 
# @param frag Matrix or Fragility object. Either a matrix with row as Electrode names and Column as fragility index, or a Fragility object from \code{calc_adj_frag}
# @param elecsoz Integer.  Vector soz electrodes (for good electrodes)
# 
# 
# @return list of 5 items with quantile matrix, mean and sdv from both electrodes groups
# @export
# 
# @examples
# data("pt01Frag")
# data("pt01Epoch")
# sozindex<-attr(pt01Epoch,"sozindex")
# fragstat<-frag_stat(frag=pt01Frag, elecsoz=sozindex)
# frag_stat <- function(frag, elecsoz){
#   if (is(frag, "Fragility")) {
#     frag <- frag$frag
#   }
#   
#   n_elec <- nrow(frag)    # electrode number
#   n_steps <- ncol(frag)   # Time window number
#   
#   elecvec=1:n_elec
#   elecsozc=which(!elecvec%in%elecsoz)
# 
#   # create separate heatmaps for soz and sozc for quantile calcs
#   hmapsoz <- frag[elecsoz,,drop=FALSE]
#   hmapsozc <- frag[elecsozc,,drop=FALSE]
# 
#   quantilematrixsozsozc <- matrix(0,20,n_steps)
#   cmeansoz <- c(1:n_steps)*0
#   cmeansozc <- c(1:n_steps)*0
#   csdsoz <- c(1:n_steps)*0
#   csdsozc <- c(1:n_steps)*0
#   
#   
#   for(i in 1:n_steps){
#     
#     colsoz <- hmapsoz[,i]
#     colsozc <- hmapsozc[,i]
#     
#     meansoz <- mean(colsoz)
#     sdsoz <- sd(colsoz)
#     meansozc <- mean(colsozc)
#     sdsozc <- sd(colsozc)
#     
#     cmeansoz[i] <- meansoz
#     cmeansozc[i] <- meansozc
#     csdsoz[i] <- sdsoz
#     csdsozc[i] <- sdsozc
#     
#     f10colsoz<-quantile(colsoz,probs=c(0.1))
#     f20colsoz<-quantile(colsoz,probs=c(0.2))
#     f30colsoz<-quantile(colsoz,probs=c(0.3))
#     f40colsoz<-quantile(colsoz,probs=c(0.4))
#     f50colsoz<-quantile(colsoz,probs=c(0.5))
#     f60colsoz<-quantile(colsoz,probs=c(0.6))
#     f70colsoz<-quantile(colsoz,probs=c(0.7))
#     f80colsoz<-quantile(colsoz,probs=c(0.8))
#     f90colsoz<-quantile(colsoz,probs=c(0.9))
#     f100colsoz<-quantile(colsoz,probs=c(1.0))
#     
#     f10colsozc<-quantile(colsozc,probs=c(0.1))
#     f20colsozc<-quantile(colsozc,probs=c(0.2))
#     f30colsozc<-quantile(colsozc,probs=c(0.3))
#     f40colsozc<-quantile(colsozc,probs=c(0.4))
#     f50colsozc<-quantile(colsozc,probs=c(0.5))
#     f60colsozc<-quantile(colsozc,probs=c(0.6))
#     f70colsozc<-quantile(colsozc,probs=c(0.7))
#     f80colsozc<-quantile(colsozc,probs=c(0.8))
#     f90colsozc<-quantile(colsozc,probs=c(0.9))
#     f100colsozc<-quantile(colsozc,probs=c(1.0))
#     
#     quantilematrixsozsozc[1,i]=f10colsoz
#     quantilematrixsozsozc[2,i]=f20colsoz
#     quantilematrixsozsozc[3,i]=f30colsoz
#     quantilematrixsozsozc[4,i]=f40colsoz
#     quantilematrixsozsozc[5,i]=f50colsoz
#     quantilematrixsozsozc[6,i]=f60colsoz
#     quantilematrixsozsozc[7,i]=f70colsoz
#     quantilematrixsozsozc[8,i]=f80colsoz
#     quantilematrixsozsozc[9,i]=f90colsoz
#     quantilematrixsozsozc[10,i]=f100colsoz
#     quantilematrixsozsozc[11,i]=f10colsozc
#     quantilematrixsozsozc[12,i]=f20colsozc
#     quantilematrixsozsozc[13,i]=f30colsozc
#     quantilematrixsozsozc[14,i]=f40colsozc
#     quantilematrixsozsozc[15,i]=f50colsozc
#     quantilematrixsozsozc[16,i]=f60colsozc
#     quantilematrixsozsozc[17,i]=f70colsozc
#     quantilematrixsozsozc[18,i]=f80colsozc
#     quantilematrixsozsozc[19,i]=f90colsozc
#     quantilematrixsozsozc[20,i]=f100colsozc
#     
#   }
#   
#   FragStat(
#     qmatrix=quantilematrixsozsozc,
#     cmeansoz=cmeansoz,
#     cmeansozc=cmeansozc,
#     csdsoz=csdsoz,
#     csdsozc=csdsozc
#       )
# }


predictRidge <- function(xt, A) {
    ## the data matrix
    if (nrow(A) == ncol(A) + 1) {
        x <- cbind(1, as.matrix(xt))
    } else {
        x <- as.matrix(xt)
    }
    x %*% A
}
