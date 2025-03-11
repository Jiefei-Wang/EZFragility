#' fit a generalized linear model to compute adjacency matrix A
#'
#' A x(t) = x(t+1)
#'
#' @param xt matrix. iEEG time series for a given window,
#' with time points as rows and electrodes names as columns
#' @param xtp1 matrix. the iEEG time serie at the next time point,
#' with time points as rows and electrodes names as columns
#' @param lambda Numeric Vector. A user supplied lambda sequence.
#' @param intercept Boolean. Should intercept(s) be fitted (default=TRUE) or set to zero (FALSE)
#'
#' @return adjacency matrix A
ridge <- function(xt, xtp1, lambda, intercept = FALSE) {
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
ridgesearchlambdadichomotomy <- function(xt, xtp1, intercept = FALSE){
  if(!identical(dim(xt),dim(xtp1)))
    stop("Unmatched dimension")
  nel <- ncol(xt)
  ## Coefficient matrix A
  ## each column is coefficients from a linear regression
  ## formula: xtp1 = xt*A + E
  lambdamin <- 0.0001
  lambdamax <- 10


  Aa <- ridge(xt,xtp1,lambda=lambdamin,intercept=F)

  stableam <- TRUE

  nel <- ncol(Aa)
  e <- Mod(eigen(Aa)$values)
  me <- max(e)

  if (me >= 1) {
    stableam <- FALSE
  }

  #print(stableam)

  if(stableam){
    lambdaopt <- lambdamin
    Afin <- Aa
  }else{

    stablea <- stableam
    lambdaa <- lambdamin

    lambdab <- lambdamax

    Ab<-ridge(xt,xtp1,lambda=lambdab,intercept=F)

    stableb <- TRUE

    nel <- ncol(Ab)
    e <- Mod(eigen(Ab)$values)
    me <- max(e)

    if (me >= 1) {
      stableb <- FALSE
    }

    #print(stableb)
    k <- 0
    while(k<20){
      lambdac <- (lambdaa + lambdab)*0.5

      Ac<-ridge(xt,xtp1,lambda=lambdac,intercept=F)

      stablec <- TRUE

      nel <- ncol(Ac)
      e <- Mod(eigen(Ac)$values)
      me <- max(e)

      if (me >= 1) {
        stablec <- FALSE
      }

      if(!stablec){
        lambdaa <- lambdac
        lambdaopt <- lambdab

      }else{
        lambdab <- lambdac
        lambdaopt <- lambdac
      }
      k <- k+1

      # print("ite")
      # print(k)
      # print(lambdac)
      # print(stablec)
      # print(lambdaopt)
    }
  }

  Afin <- ridge(xt,xtp1,lambda=lambdaopt,intercept=F)

  attr(Afin, "lambdaopt") <- lambdaopt
  Afin
}

