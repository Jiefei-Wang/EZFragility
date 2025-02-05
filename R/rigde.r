RIDGE <- function(x, x1, l) {
  n <- ncol(x)
  cb <- \(i) { glmnet::glmnet(x = x, y = x1[, i], lambda = l,
    nlambda = 1, alpha = 0, standardize = FALSE, intercept = FALSE,
    family = "gaussian", type.gaussian = "covariance" )[["beta"]]@x
  }
  A <- vapply(seq_len(n), cb, numeric(n), USE.NAMES = FALSE)
  isStable <- (eigen(A, FALSE, TRUE)$values |> Mod() |> max()) < 1.0
  structure(A, lambda = l, stable = isStable)
}

RIDGESearchLambdaDichomotomy <- function(xt, xtp1) {
  if (!identical(dim(xt), dim(xtp1))) stop("Unmatched dimension")
  e <- environment()
  low = lmbd <- 1e-4
  high <- 10.0
  A <- RIDGE(xt, xtp1, lmbd)
  if (!attr(A, "stable")) {
    for (i in seq_len(20L)) {
      lmbd <- (low + high) * .5
      A <- RIDGE(xt, xtp1, lmbd)
      if (attr(A, "stable")) high <- lmbd
      else low <- lmbd
      if (high - low < 0.01) break
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
  if (!identical(dim(xt), dim(xtp1))) {
    stop("Unmatched dimension")
  }
  nel <- ncol(xt)
  ## Coefficient matrix A
  ## each column is coefficients from a linear regression
  ## formula: xtp1 = xt*A + E
  A <- matrix(0, nel + intercept, nel)
  ## for each electrode
  for (i in seq_len(nel)) {
    y <- xtp1[, i]
    fit <- glmnet::glmnet(xt, y,
                  alpha = 0, lambda = lambda,
                  standardize = FALSE, intercept = intercept
    )
    # fit <- glmnet::cv.glmnet(xt, y,
    #               alpha = 0,
    #               standardize = FALSE, intercept = intercept
    # )

    if (intercept) {
      A[, i] <- as.numeric(coef(fit))
    } else {
      A[, i] <- coef(fit)[-1]
    }
  }
  AEigen <- eigen(A)
  e <- Mod(AEigen$values)

  A
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

