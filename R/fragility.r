#' Compute the normalized fragility row for adjacency matrix A
#'
#' @param A 
#' @param nSearch 
#'

fragilityRowNormalized <- function(A, nSearch = 100) {
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
  omvec <- seq(0, 1, length.out = nSearch + 1)[-1]

  b <- c(0, -1)
  ## for each electrode
  for (i in 1:nel) {
    ## indicate which electrode is disturbed (ith)
    ek <- rep(0, nel)
    ek[i] <- 1
    tek <- t(ek)
    minNorm <- 100000
    minPerturbColumn <- NA
    for (k in seq_len(nSearch)) {
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

#' Title
#'
#' @param A 
#' @param nSearch 
#'
fragilityRow <- function(A, nSearch = 100) {
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
  omvec <- seq(0, 1, length.out = nSearch + 1)[-1]

  b <- c(0, -1)
  ## for each electrode
  for (i in 1:nel) {
    ## indicate which electrode is disturbed (ith)
    ek <- rep(0, nel)
    ek[i] <- 1
    tek <- t(ek)
    minNorm <- 100000
    minPerturbColumn <- NA
    for (k in seq_len(nSearch)) {
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


frag_quantile <- function(repository, f, t_window, t_step, soz, sozc){
  n_tps <- length(repository$voltage$dimnames$Time)
  n_elec <- length(repository$voltage$dimnames$Electrode)
  n_steps <- floor((n_tps - t_window) / t_step) + 1
  epoch_time_window <- repository$time_windows[[1]]
  fs <- round(repository$sample_rate,-1)
  if(any(repository$electrode_table$Label == "NoLabel")) {
    elec_names <- repository$electrode_table$Electrode[match(c(soz,sozc), repository$electrode_table$Electrode)]
    elec_names <- as.character(elec_names)
  } else {
    elec_names <- repository$electrode_table$Label[match(c(soz,sozc), repository$electrode_table$Electrode)]
  }

  # create fragility map with soz electrodes separated from sozc electrodes
  fmap <- f[as.character(c(soz,sozc)),]
  stimes <- (seq_len(n_steps)-1)*t_step/fs+epoch_time_window[1]

  # raw fragility map
  fplot <- expand.grid(Time = stimes, Electrode = elec_names)
  fplot$Value <- c(t(fmap))

  # create separate heatmaps for soz and sozc for quantile calcs
  hmapsoz <- fmap[as.character(soz),]
  hmapsozc <- fmap[as.character(sozc),]

  #f90soz=quantile(hmapsoz, probs=c(0.9))
  #f90sozc=quantile(hmapsozc,probs=c(0.9))
  #interpretabilityratiosoz=f90soz/f90sozc

  quantilematrixsozsozc=matrix(0,20,length(stimes))
  cmeansoz=c(1:length(stimes))*0
  cmeansozc=c(1:length(stimes))*0
  csdsoz=c(1:length(stimes))*0
  csdsozc=c(1:length(stimes))*0

  for(i in 1:length(stimes)){

    colsoz=hmapsoz[,i]
    colsozc=hmapsozc[,i]

    meansoz=mean(colsoz)
    sdsoz=sd(colsoz)
    meansozc=mean(colsozc)
    sdsozc=sd(colsozc)

    cmeansoz[i]=meansoz
    cmeansozc[i]=meansozc
    csdsoz[i]=sdsoz
    csdsozc[i]=sdsozc

    f10colsoz<-quantile(colsoz,probs=c(0.1))
    f20colsoz<-quantile(colsoz,probs=c(0.2))
    f30colsoz<-quantile(colsoz,probs=c(0.3))
    f40colsoz<-quantile(colsoz,probs=c(0.4))
    f50colsoz<-quantile(colsoz,probs=c(0.5))
    f60colsoz<-quantile(colsoz,probs=c(0.6))
    f70colsoz<-quantile(colsoz,probs=c(0.7))
    f80colsoz<-quantile(colsoz,probs=c(0.8))
    f90colsoz<-quantile(colsoz,probs=c(0.9))
    f100colsoz<-quantile(colsoz,probs=c(1.0))

    f10colsozc<-quantile(colsozc,probs=c(0.1))
    f20colsozc<-quantile(colsozc,probs=c(0.2))
    f30colsozc<-quantile(colsozc,probs=c(0.3))
    f40colsozc<-quantile(colsozc,probs=c(0.4))
    f50colsozc<-quantile(colsozc,probs=c(0.5))
    f60colsozc<-quantile(colsozc,probs=c(0.6))
    f70colsozc<-quantile(colsozc,probs=c(0.7))
    f80colsozc<-quantile(colsozc,probs=c(0.8))
    f90colsozc<-quantile(colsozc,probs=c(0.9))
    f100colsozc<-quantile(colsozc,probs=c(1.0))

    quantilematrixsozsozc[1,i]=f10colsoz
    quantilematrixsozsozc[2,i]=f20colsoz
    quantilematrixsozsozc[3,i]=f30colsoz
    quantilematrixsozsozc[4,i]=f40colsoz
    quantilematrixsozsozc[5,i]=f50colsoz
    quantilematrixsozsozc[6,i]=f60colsoz
    quantilematrixsozsozc[7,i]=f70colsoz
    quantilematrixsozsozc[8,i]=f80colsoz
    quantilematrixsozsozc[9,i]=f90colsoz
    quantilematrixsozsozc[10,i]=f100colsoz
    quantilematrixsozsozc[11,i]=f10colsozc
    quantilematrixsozsozc[12,i]=f20colsozc
    quantilematrixsozsozc[13,i]=f30colsozc
    quantilematrixsozsozc[14,i]=f40colsozc
    quantilematrixsozsozc[15,i]=f50colsozc
    quantilematrixsozsozc[16,i]=f60colsozc
    quantilematrixsozsozc[17,i]=f70colsozc
    quantilematrixsozsozc[18,i]=f80colsozc
    quantilematrixsozsozc[19,i]=f90colsozc
    quantilematrixsozsozc[20,i]=f100colsozc

  }

  quantilesname<-c("SOZ(10th)","SOZ(20th)","SOZ(30th)","SOZ(40th)","SOZ(50th)",
                   "SOZ(60th)","SOZ(70th)","SOZ(80th)","SOZ(90th)","SOZ(100th)",
                   "SOZc(10th)","SOZc(20th)","SOZc(30th)","SOZc(40th)","SOZc(50th)",
                   "SOZc(60th)","SOZc(70th)","SOZc(80th)","SOZc(90th)","SOZc(100th)")
  quantileplot<- expand.grid(Time = stimes, Stats=quantilesname)
  quantileplot$Value <- c(t(quantilematrixsozsozc))

  dimnames(quantilematrixsozsozc) <- list(
    Quantile = quantilesname,
    Time = stimes
  )

  return(list(
    fplot = fplot,
    q_matrix = quantilematrixsozsozc,
    q_plot = quantileplot
  ))
}

mean_f_calc <- function(repository, f, soz, sozc) {
  mean_f_soz <- rep(0,dim(f)[2])
  mean_f_sozc <- rep(0,dim(f)[2])
  se_f_soz <- rep(0,dim(f)[2])
  se_f_sozc <- rep(0,dim(f)[2])
  for (i in seq_len(dim(f)[2])){
    mean_f_soz[i] <- mean(f[as.character(soz),i])
    se_f_soz[i] <- sd(f[as.character(soz),i])/sqrt(length(soz))
    mean_f_sozc[i] <- mean(f[as.character(sozc),i])
    se_f_sozc[i] <- sd(f[as.character(sozc),i]/sqrt(length(sozc)))
  }
  return(list(
    mean_f_soz = mean_f_soz,
    mean_f_sozc = mean_f_sozc,
    se_f_soz = se_f_soz,
    se_f_sozc = se_f_sozc
  ))
}

threshold_buckets <- function(mat, thresholds) {

  thresholds <- sort(thresholds)

  if (max(mat) > 1) {
    stop("Matrix values must be between 0 and 1!")
  }

  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      k <- 1
      while(mat[i,j] >= thresholds[k]) {
        k <- k + 1
      }
      mat[i,j] <- mean(thresholds[(k-1):k])
    }
  }

  mat <- (mat - min(mat))/(max(mat) - min(mat)) # normalize to between 0 and 1

  return(mat)
}

moving_average <- function(x, n) {
  x_ma <- stats::filter(x, rep(1 / n, n), sides = 2)
  x_ma[is.na(x_ma)] <- x[is.na(x_ma)] # replace NAs at beginning and end with original values
  return(as.vector(x_ma))
}

threshold_fragility <- function(repository, adj_frag_info, t_step, threshold_start, threshold_end, threshold = 0.5) {
  n_windows <- dim(adj_frag_info$adj)[3]

  # convert from input t_start and t_end to timewindow indices
  tw_start <- floor(which.min(abs(threshold_start-repository$voltage$dimnames$Time))/t_step)
  tw_end <- floor(which.min(abs(threshold_end-repository$voltage$dimnames$Time))/t_step)
  if (tw_end > n_windows) { tw_end <- n_windows }

  # subset fragility matrix to specified timewindows
  mat <- adj_frag_info$frag[,tw_start:tw_end]

  avg_f <- rowMeans(mat)
  elec <- which(avg_f > threshold)

  return(list(
    avg_f = avg_f,
    elecnames = attr(elec, "names")
  ))
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