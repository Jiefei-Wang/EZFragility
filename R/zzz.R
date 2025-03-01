#' @import methods
#' @import stats
NULL

pkgData <- new.env()
pkgData$debug <- FALSE

debug <- function(){
    pkgData$debug <- TRUE
}

undebug <- function(){
    pkgData$debug <- FALSE
}