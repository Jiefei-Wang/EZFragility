isWholeNumber <- function(x) {
  return(x %% 1 == 0)
}

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

