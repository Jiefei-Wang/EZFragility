#' @import methods
#' @import stats
#' @import Epoch
#' @importFrom rlang .data
#' @importFrom glue glue
#' @importFrom foreach foreach %dopar%
#' @importFrom progress progress_bar
#' @importFrom ramify pprint
#' @importFrom ggplot2 ggplot aes geom_line xlab ylab geom_ribbon scale_color_manual theme element_text scale_y_continuous theme_minimal scale_x_discrete geom_tile labs
#' @importFrom ggtext element_markdown
#' @importFrom viridis scale_fill_viridis
NULL

pkgData <- new.env()
pkgData$debug <- FALSE

debug <- function() {
    pkgData$debug <- TRUE
}

undebug <- function() {
    pkgData$debug <- FALSE
}
