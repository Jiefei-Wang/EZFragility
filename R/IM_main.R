# Replaces calc_adj_frag (I don't like snake_case - I prefer CamelCase)
GetAdjFrag <- function(ieegts, t_window, t_step, lambda = NULL, nSearch = 10L) {
  stopifnot(isWholeNumber(t_window))
  stopifnot(isWholeNumber(t_step))
  stopifnot(is.null(lambda) | is.numeric(lambda))
  stopifnot(nrow(ieegts) >= t_window)
  self <- environment()
  ConLogger <- ConsoleLogger(self)
  scaling <- 10^floor(log10(max(ieegts)))
  ieegts  <- ieegts / scaling
  # Electrode count and names
  elCnt <- ncol(ieegts)
  elNms <- colnames(ieegts)
  # Number/sequence of steps
  nsteps  <- floor((nrow(ieegts) - t_window) / t_step) + 1L
  STEPS   <- seq_len(nsteps)
  # Pre-allocate output
  dm   <- c(elCnt, elCnt, nsteps)
  dmn  <- list(Electrode  = elNms, Step = STEPS)
  dmnA <- list(Electrode1 = elNms, Electrode2 = elNms, Step = STEPS)
  A    <- array(.0, dim = dm,     dimnames = dmnA)
  R2   <- array(.0, dim = dm[-1], dimnames = dmn)
  f = fR <- R2
  lbd <- rep(0, nsteps) |> setNames(STEPS)
  # Indices of window at time 0
  i0 <- seq_len(t_window - 1L)
  ConLogger$ProcessInfo()
  for (iw in STEPS) {
    si   <- i0 + (iw - 1L) * t_step
    xt   <- ieegts[si, ]
    xtp1 <- ieegts[si + 1L, ]
    ConLogger$RidgeStart()
    adjMatrix <- RIDGESearchLambdaDichomotomy(xt, xtp1)
    A[,, iw]  <- adjMatrix
    R2[, iw]  <- RidgeRSQ(xt, xtp1, adjMatrix)
    ConLogger$FragStart()
    f[,  iw]  <- FragilityROW(adjMatrix, nSearch)
    fR[, iw]  <- rank(f[, iw]) / elCnt # ranks should probably be here...
    lbd[[iw]] <- attr(adjMatrix, "lambda")
    ConLogger$StepEnd()
  }
  ConLogger$TotalTime()
  Fragility(
    ieegts = ieegts,
    adj = A,
    R2 = R2,
    frag = f,
    frag_ranked = fR,
    lambdas = lbd
  )
}

# Optional utility which prints info while the above is running.
# I wrote it because I was bored to death waiting for the thing to run
# Also, it helped me debug code bottlenecks.
ConsoleLogger <- \(e) {
  start <- Sys.time()
  stepStart = TimeKeeper = RidgeRun = FragRun <-  NULL;
  InfoDash = TabDash <- NULL
  t_window <- e$t_window;
  t_step <- e$t_step;
  samples <- nrow(e$ieegts)
  self <- environment()
  # Shifts to the right all strings of a list with a number of blanks
  shift <- \(strL, nBlanks) {
    pre <- paste(rep(" ", nBlanks), collapse = "")
    lapply(strL, \(x) sprintf("%s%s", pre, x)) |> unlist()
  }
  # Get number of seconds from a reference timestamp in specified format
  getTimeSecs <- \(ref) {
    difT <- difftime(Sys.time(), ref, units = "secs") |> as.double()
    sprintf("%.2f", difT)
  }
  # Save time of Ridge (and step) start and print step number
  RidgeStart <- \() {
    self$stepStart <- Sys.time()
    self$TimeKeeper <- self$stepStart
    sprintf("%7d ", e$iw) |> cat()
  }
  # Save time of Fragility start and print Ridge runtime
  FragStart <- \() {
    self$RidgeRun <- getTimeSecs(self$TimeKeeper)
    self$TimeKeeper <- Sys.time()
    sprintf("%11s ", self$RidgeRun) |> cat()
  }
  # Print runtime for Fragility and for the whole step
  StepEnd <- \() {
    sprintf("%11s ", getTimeSecs(self$TimeKeeper)) |> cat()
    sprintf("%7s",   getTimeSecs(self$stepStart))  |> cat("\n")
  }
  # Print the total runtime of the whole process
  TotalTime <- \() {
    total <- Sys.time() - start
    fmt <- "  Total Runtime:  %21.2f %s"
    self$StepDash |> shift(2) |> cat("\n")
    sprintf(fmt, as.double(total), attr(total, "units")) |> cat("\n")
  }
  # Prints process specifications
  initTab <- \(nb = 2) {
    header <- list("Samples", "Window", "Shift", "Steps")
    values <- list(samples, t_window, t_step, e$nsteps)
    hFmt <- do.call(sprintf, c(list(" %7s |%7s |%6s |%6s "), header))
    vFmt <- do.call(sprintf, c(list(" %7d  %7d  %6d  %6d "), values))
    self$InfoDash <- paste(rep("-", nchar(hFmt)), collapse = "")
    c(self$InfoDash, hFmt, vFmt) |> shift(nb) |> cat(sep = "\n")
  }
  # Prints the headers of the step time table
  stepTab <- \(nb = 2) {
    header <- c("Step", "Adjacency", "Fragility", "Total")
    hFmt <- do.call(sprintf, c(list("%5s | %9s | %9s | %5s"), header))
    len <- nchar(hFmt)
    self$StepDash <- paste(rep("-", len + 1), collapse = "")
    subTitle <- paste(rep("-", len - 12), collapse = "") 
    subtext <- "Runtime (sec)"
    id <- (1 + 0.5 * (nchar(subTitle) - nchar(subtext))) |> ceiling()
    substr(subTitle, id, id + nchar(subtext)) <- subtext
    c(self$StepDash, shift(subTitle, 6), hFmt) |> shift(nb) |> cat(sep = "\n")
  }
  ProcessInfo <- \() { initTab(); stepTab() }
  self
}

