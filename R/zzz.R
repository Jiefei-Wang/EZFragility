#' @import methods
#' @import stats
#' @importFrom rlang .data
#' @importFrom glue glue
#' @importFrom foreach foreach %dopar%
#' @importFrom progress progress_bar
#' @importFrom ramify pprint
NULL

pkgData <- new.env()
pkgData$debug <- FALSE

debug <- function() {
    pkgData$debug <- TRUE
}

undebug <- function() {
    pkgData$debug <- FALSE
}
