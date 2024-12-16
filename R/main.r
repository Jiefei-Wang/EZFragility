#' Calculate adjacency matrices and fragility heatmap from iEEG recording
#' 
#' A very long long description to give the details of the function
#' 
#' @param readings Numeric. A matrix of voltage readings, 
#' with time points as rows and electrodes as columns
#' @param t_window Integer. The number of time points to use in each window
#' @param t_step Integer.The number of time points to move the window each time
#' @param lambda Numeric. The lambda value to use in the ridge regression. 
#' If NULL, the lambda will be chosen automatically
#' 
#' @return A list containing the normalized readings, 
#' adjacency matrices, fragility, and R^2 values
#' 
#' @examples
#' data <- matrix(rnorm(100), nrow = 10)
#' t_window <- 10
#' t_step <- 5
#' lambda <- 0.1
#' calc_adj_frag(readings = data, t_window = t_window, t_step = t_step, lambda = lambda)
#' 
#' @export 
calc_adj_frag <- function(readings, t_window, t_step, lambda = NULL) {
    ## check the input types
    stopifnot(isWholeNumber(t_window))
    stopifnot(isWholeNumber(t_step))
    stopifnot(is.null(lambda) | is.numeric(lambda))

    ## The input matrix must have at least t_window rows
    stopifnot(nrow(readings) >= t_window)


    ## Number of electrodes and time points
    n_tps <- nrow(readings)
    n_elec <- ncol(readings)

    electrode_list <- colnames(readings)

    # Number of steps
    n_steps <- floor((n_tps - t_window) / t_step) + 1

    scaling <- 10^floor(log10(max(readings)))
    readings <- readings / scaling

    ## create adjacency array (array of adj matrices for each time window)
    ## iw: The index of the window we are going to calculate fragility
    res <- lapply(seq_len(n_steps), function(iw) {
        ## Sample indices for the selected window
        si <- seq_len(t_window - 1) + (iw - 1) * t_step
        ## measurements at time point t
        xt <- readings[si, ]
        ## measurements at time point t plus 1
        xtp1 <- readings[si + 1, ]

        ## Coefficient matrix A (adjacency matrix)
        ## each column is coefficients from a linear regression
        ## formula: xtp1 = xt*A + E
        if (is.null(lambda)) {
            Ai <- ridgesearchlambdadichomotomy(xt, xtp1, intercept = FALSE, iw = iw)
        } else {
            Ai <- ridge(xt, xtp1, intercept = FALSE, lambda = lambda, iw = iw)
        }

        R2 <- ridgeR2(xt, xtp1, Ai)

        list(Ai = Ai, R2 = R2)
    })

    A <- unlist(lapply(res, function(w) {
        w$Ai
    }))
    ## TODO: Why do you want to do this? very error prone
    dim(A) <- c(n_elec, n_elec, n_steps)
    dimnames(A) <- list(
        Electrode1 = electrode_list,
        Electrode2 = electrode_list,
        Step = seq_len(n_steps)
    )

    R2 <- unlist(lapply(res, function(w) {
        w$R2
    }))
    dim(R2) <- c(n_elec, n_steps)
    dimnames(R2) <- list(
        Electrode = electrode_list,
        Step = seq_len(n_steps)
    )

    if (is.null(lambda)){
        lambdas <- sapply(res, function(w) {
            attr(w$Ai, "lambdaopt")
        })
    } else {
        lambdas <- rep(lambda, length(res))
    }
    

    # calculate fragility
    f <- sapply(seq_len(n_steps), function(iw) {
        fragilityRowNormalized(A[, , iw])
    })
    dimnames(f) <- list(
        Electrode = electrode_list,
        Step = seq_len(n_steps)
    )

    ## TODO: Is this consistent with the method in the paper?
    # ranked fragility map
    f_rank <- matrix(rank(f), nrow(f), ncol(f))
    attributes(f_rank) <- attributes(f)
    f_rank <- f_rank / max(f_rank)

    return(list(
        voltage = readings, 
        adj = A,
        frag = f,
        frag_ranked = f_rank,
        R2 = R2,
        lambdas = lambdas
    ))
}
