#' Visualization functions (raw signal, fragility matrix)
#'
#' plot fragility heatmaps with electrodes marked as soz colored
#'
#' @param frag Fragility object from \code{calcAdjFrag}
#' @param sozIndex Integer or string. A group of electrodes to mark as in the Seizure Onset Zone (SOZ)
#' @param title String. Figure title
#' @param display Integer or string. Vector electrodes to display
#'
#' @return Heatmap plot of the fragility matrix with soz electrodes in blue in the bottom
#'
#' @examples
#' # use integer index for display and soz electrodes
#' data("pt01Epoch")
#' sozIndex <- attr(pt01Epoch, "sozIndex")
#' data("pt01Frag")
#' timeRange <- c(-1, 2)
#' display <- c(sozIndex, 77:80)
#' fragplot <- heatmapFrag(
#'     frag = pt01Frag, sozID = sozIndex,
#'     timeRange = timeRange, title = "PT01 seizure 1", display = display
#' )
#' fragplot
#'
#'
#' # use electrodes name for display and soz electrodes
#' data("pt01Epoch")
#' sozNames <- attr(pt01Epoch, "sozNames")
#' data("pt01Frag")
#' timeRange <- c(-1, 2)
#' display <- c(sozNames, "MLT1", "MLT2", "MLT3", "MLT4")
#' fragplot <- heatmapFrag(
#'     frag = pt01Frag, sozID = sozNames,
#'     timeRange = timeRange, title = "PT01 seizure 1", display = display
#' )
#' fragplot
#'
#' # save plot to file with ggplot2
#' data("pt01Epoch")
#' data("pt01Frag")
#' sozIndex <- attr(pt01Epoch, "sozIndex")
#' timeRange <- c(-10, 10)
#' display <- c(sozIndex, 77:80)
#' pathplot <- "~"
#' title <- "PT01sz1"
#' resfile <- paste(pathplot, "/FragilityHeatMap", title, ".png", sep = "")
#' fragplot <- heatmapFrag(
#'     frag = pt01Frag, sozID = sozIndex, timeRange = timeRange,
#'     title = title, display = display
#' )
#' fragplot
#' ggplot2::ggsave(resfile)
#'
#' @export
heatmapFrag <- function(
    frag, 
    sozIndex = NULL) {
    ## TODO: make sozID an optional
    ## TODO: add plot support to frag
    fragMat <- frag$frag

    elecNames <- rownames(fragMat)
    sozIndex <- checkIndex(sozIndex, elecNames)

    nElec <- nrow(fragMat)
    nTimes <- ncol(fragMat)
    elecNames <- rownames(fragMat)
    group1 <- sozIndex
    group2 <- setdiff(seq_len(nElec), sozIndex)

    elecColor <- rep("black", nElec)
    elecColor[seq_along(group1)] <- "blue"

    startTime <- frag$startTimes
    if (is.null(startTime)) {
        xlabel <- "Time Index"
        stimes <- seq_len(nTimes)
    } else {
        xlabel <- "Time (s)"
        stimes <- startTime
    }

    ## prepare the data.frame for visualization
    allIndex <- c(group1, group2)
    fragMatReorderd <- fragMat[allIndex, ]

    rownames(fragMatReorderd) <- frag$electrodes[allIndex]
    colnames(fragMatReorderd) <- stimes

    ## make fragDf a long format data.frame
    ## three columns: Electrode, Time, Value
    df_long <- as.data.frame(as.table(fragMatReorderd))
    colnames(df_long) <- c("Electrode", "Time", "Value")

    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data$Time, y = .data$Electrode, fill = .data$Value)) +
        ggplot2::geom_tile() +
        ggplot2::theme(plot.title = ggtext::element_markdown(hjust = 0.5)) +
        ggplot2::labs(x = xlabel, y = "Electrode", size = 2) +
        viridis::scale_fill_viridis(option = "turbo") +
        ggplot2::geom_vline(
            xintercept = 0,
            color = "black", linetype = "dashed", linewidth = 1
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggtext::element_markdown(size = 6, colour = elecColor), # Adjust depending on electrodes
        )

    return(p)
}

#' Visualization of ictal iEEG
#'
#' @inheritParams heatmapFrag
#' @param ieegts Matrix or Fragility object. Either a matrix of iEEG time
#' series x(t), with time points as rows and electrodes names as columns,
#' or a Fragility object from \code{calcAdjFrag}
#' @param title String. Figure title
#' @param display Integer or string. Vector electrodes to display
#' @return plot raw signal
#'
#' @examples
#' data("pt01Epoch")
#' sozIndex <- attr(pt01Epoch, "sozIndex")
#' display <- c(sozIndex, 77:80)
#' timeRange <- c(-1, 2)
#' iEEGplot <- visuIEEGData(ieegts = pt01Epoch, timeRange = timeRange, display = display)
#' iEEGplot
#' @export
visuIEEGData <- function(ieegts, timeRange = NULL, title = "Patient name seizure number", display = NULL) {
    titlepng <- title
    if (is.null(display)) {
        display <- seq_len(ncol(ieegts))
    }

    elecNames <- colnames(ieegts)
    displayid <- checkIndex(display, elecNames)

    scaling <- 10^floor(log10(max(ieegts)))
    plotData <- ieegts[, displayid] / scaling
    gaps <- 2
    displayNames <- colnames(ieegts)[displayid]
    nElec <- length(displayid)
    nt <- nrow(plotData)
    if (is.null(timeRange)) {
        xlabel <- "Time Index"
        stimes <- seq_len(nt)
    } else {
        xlabel <- "Time (s)"
        stimes <- seq(timeRange[1], timeRange[2], length.out = nt)
    }


    for (i in 1:ncol(plotData)) {
        plotData[, i] <- (plotData[, i] - mean(plotData[, i])) +
            (ncol(plotData) - i) * gaps
    }
    plotData <- as.data.frame(plotData)
    plotData$stimes <- stimes
    breakplot <- (c(1:nElec) - 1) * gaps

    p <- ggplot2::ggplot(data = plotData, ggplot2::aes(x = .data$stimes, y = .data$plotData)) +
        ggplot2::ggtitle(titlepng) +
        ggplot2::labs(x = xlabel, y = "Electrode", size = 2) +
        ggplot2::geom_vline(
            xintercept = 0,
            color = "black", linetype = "dashed", linewidth = 1
        )

    for (i in displayNames) {
        p <- p + ggplot2::geom_line(ggplot2::aes(y = .data[[i]]))
    }
    displayNames <- rev(displayNames)
    p <- p + ggplot2::scale_y_continuous(labels = displayNames, breaks = breakplot)

    return(p)
}

#' Plot Fragility time quantiles for two electrodes group marked as soz non marked as soz
#'
#' @inheritParams heatmapFrag
#' @param FragStatObj Matrix or FragStat object, either a quantile matrix
#' for the two groups or a FragStat object from \code{fragStat}
#' @param title String. Figure title
#'
#' @return Quantile plot
#' @export
#'
#' @examples
#'
#' timeRange <- c(-1, 2)
#' data("pt01Epoch")
#' sozIndex <- attr(pt01Epoch, "sozIndex")
#' data("pt01Frag")
#' # compute fragility statistics evolution with time (mean and standard deviation) for soz and
#' # non soz groups
#' pt01fragstat <- fragStat(frag = pt01Frag, sozID = sozIndex)
#' plotFragQuantile(FragStatObj = pt01fragstat, timeRange = timeRange)
plotFragQuantile <- function(FragStatObj, timeRange = NULL, title = "Fragility Quantiles over time") {
    if (is(FragStatObj, "FragStat")) {
        qmatrix <- FragStatObj$qmatrix
    }

    nw <- ncol(qmatrix)
    if (is.null(timeRange)) {
        xlabel <- "Time Index"
        stimes <- seq_len(nw)
    } else {
        xlabel <- "Time (s)"
        stimes <- seq(timeRange[1], timeRange[2], length.out = nw)
    }


    quantilesName <- rownames(qmatrix)
    quantilePlot <- expand.grid(Time = stimes, Stats = quantilesName)
    quantilePlot$Value <- c(t(qmatrix))

    titlepng <- title

    p <- ggplot2::ggplot(quantilePlot, ggplot2::aes(x = .data$Time, y = .data$Stats, fill = .data$Value)) +
        ggplot2::geom_tile() +
        ggplot2::ggtitle(titlepng) +
        ggplot2::labs(x = xlabel, y = "Quantiles", size = 2) +
        viridis::scale_fill_viridis(option = "turbo") + #

        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 4), # Adjust depending on electrodes
        )

    if (!is.null(timeRange)) {
        p <- p + ggplot2::geom_vline(
            xintercept = 0,
            color = "black", linetype = "dashed", linewidth = 1
        )
    }

    return(p)
}


#'  Plot Fragility time distribution for two electrodes group marked and non-marked as soz
#'
#' @inheritParams heatmapFrag
#' @param stat FragStat object, a FragStat object from \code{fragStat}
#' @param title String. Figure title
#'
#' @return plot fragility distribution
#' @export
#'
#' @examples
#' data("pt01Epoch")
#' sozindex <- attr(pt01Epoch, "sozindex")
#' # Load the precomputed fragility object
#' timeRange <- c(-10, 10)
#' data("pt01Frag")
#' # compute fragility statistics evolution with time (mean and standard deviation) for soz and
#' # non soz groups
#' pt01fragstat <- fragStat(frag = pt01Frag, sozID = sozindex)
#' # plot the statistical results
#' pfragstat <- plotFragDistribution(stat = pt01fragstat, timeRange = timeRange)
#' pfragstat
plotFragDistribution <- function(
    stat = NULL,
    timeRange = NULL,
    title = "Average Fragility over time") {
    stopifnot(is(stat, "FragStat"))
    cmeansoz <- stat$cmeansoz
    cmeansozc <- stat$cmeansozc
    csdsoz <- stat$csdsoz
    csdsozc <- stat$csdsozc

    nw <- length(cmeansoz)
    if (is.null(timeRange)) {
        xlabel <- "Time Index"
        stimes <- seq_len(nw)
    } else {
        xlabel <- "Time (s)"
        stimes <- seq(timeRange[1], timeRange[2], length.out = nw)
    }


    sozsdp <- cmeansoz + csdsoz
    sozsdm <- cmeansoz - csdsoz
    sozcsdp <- cmeansozc + csdsozc
    sozcsdm <- cmeansozc - csdsozc

    plotmeanstd <- as.data.frame(stimes)
    colnames(plotmeanstd) <- "times"
    plotmeanstd$meansoz <- cmeansoz
    plotmeanstd$sozsdp <- sozsdp
    plotmeanstd$sozsdm <- sozsdm
    plotmeanstd$meansozc <- cmeansozc
    plotmeanstd$sozcsdp <- sozcsdp
    plotmeanstd$sozcsdm <- sozcsdm

    titlepng <- title
    colors <- c("SOZ +/- sd" = "red", "SOZc +/- sd" = "black")
    ggplot2::theme_grey(base_size = 22)
    p <- ggplot2::ggplot(plotmeanstd, ggplot2::aes(x = .data$times, y = .data$meansoz)) +
        ggplot2::xlab(xlabel) +
        ggplot2::ylab("Fragility") +
        ggplot2::ggtitle(titlepng) +
        ggplot2::geom_line(ggplot2::aes(y = .data$meansoz, color = "SOZ +/- sd")) +
        ggplot2::geom_line(ggplot2::aes(y = .data$sozsdp), color = "red", linetype = "dotted") +
        ggplot2::geom_line(ggplot2::aes(y = .data$sozsdm), color = "red", linetype = "dotted") +
        ggplot2::geom_line(ggplot2::aes(y = .data$meansozc, color = "SOZc +/- sd")) +
        ggplot2::geom_line(ggplot2::aes(y = .data$sozcsdp), color = "black", linetype = "dotted") +
        ggplot2::geom_line(ggplot2::aes(y = .data$sozcsdm), color = "black", linetype = "dotted") +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$sozsdm, ymax = .data$sozsdp), fill = "red", alpha = 0.5) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$sozcsdm, ymax = .data$sozcsdp), fill = "black", alpha = 0.5) +
        ggplot2::scale_color_manual(name = "Electrode groups", values = c(colors))

    ## add vertical line at time 0 if timeRange is specified
    if (!is.null(timeRange)) {
        p <- p + ggplot2::geom_vline(
            xintercept = 0,
            color = "black", linetype = "dashed", linewidth = 1
        )
    }
    return(p)
}
