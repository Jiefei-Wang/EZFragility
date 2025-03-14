#' @import methods
#' @import stats
#' @importFrom rlang .data
#' @importFrom glue glue
#' @importFrom foreach foreach %dopar%
NULL

pkgData <- new.env()
pkgData$debug <- FALSE

debug <- function(){
    pkgData$debug <- TRUE
}

undebug <- function(){
    pkgData$debug <- FALSE
}