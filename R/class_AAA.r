## This file will take precedence over all the other class files
setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("arrayOrNULL", c("array", "NULL"))
setClassUnion("numericOrNULL", c("numeric", "NULL"))

printSlotValue <- function(object, slotName, k = 3) {
  val <- methods::slot(object, slotName)
  if (is.null(val)) {
    msg <- glue::glue("{slotName}: NULL")
  } else if (is.matrix(val) || is.array(val)) {
    truncated <- if (length(as.vector(val)) > k) "..." else ""
    msg <- glue::glue(
      "{slotName} ({paste(dim(val), collapse=' x ')}): ",
      "{paste(head(as.vector(val), k), collapse=', ')}{truncated}"
    )
  } else if (is.numeric(val)) {
    if(length(val)>k){
      msg <- glue::glue(
        "{slotName} ({length(val)}): ",
        "{paste(head(val, k), collapse=', ')}..."
      )
    }else{
      msg <- 
        glue::glue(
          "{slotName}: ",
          "{paste(val, collapse=', ')}"
        )
    }
  } else {
    msg <- glue::glue("{slotName}: {val}")
  }
  cat(msg)
  cat("\n")
}





jn <- \(..., dlm = "", s = "") {
  fn <- \(x) if (length(x) > 1) paste0(x, collapse = dlm) else x
  ARGS <- lapply(list(...), fn)
  do.call(paste, c(ARGS, sep = s))
}

slotSpecs <- \(x, k = 3, dm = dim(x), vec = is.null(dm), len = length(x)) {
  val = x[seq_len(min(k, len))]
  if (is.null(x)) return(c(d = "[NULL]:", v = " NULL"))
  if (is.double(x)) val <- sprintf("%7.4f", val)
  c(d = jn("[", if (vec) len else dm, "]:", dlm = 'x'),
    v = jn(val, if (len > k) "...", dlm = ", ")
  )
}

printSlots <- \(object, nb = 0) {
  nn <- methods::slotNames(object)
  colN <- c(" Slots", " DIM|LEN", " Values")
  maxL <- sapply(colN, nchar)
  ftb <- \(i) paste0("%", -i, "s")
  meta <- list()
  dL <- nchar("DIM|LEN") 
  for (n in nn) {
    x <- c(name = n, slotSpecs(methods::slot(object, n)))
    for (i in seq_along(maxL)) maxL[i] <- max(nchar(x[[i]]), maxL[[i]])
    meta[[n]][names(x)] <- x
    dL <- max(dL, nchar(x[['d']]))
    nL <- max(nL, nchar(n))
  }
  print(maxL)
  fmts   <- c(ftb(nL), ftb(dL), "%s")
  fmt    <- list(fmt = jn(fmts, dlm = " "))
  header <- do.call(sprintf, c(fmt, " Slots", " DIM|LEN", " Values"))
  dash   <- rep("-", nchar(header)) |> jn() |> shift(nb)
  cat(dash, shift(header, nb), dash, sep = "\n")
  for (x in meta) do.call(sprintf, c(fmt, x)) |> shift(nb + 1) |> cat("\n")
  dash |> cat("\n")
}
