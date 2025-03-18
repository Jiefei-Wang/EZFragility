# A plot function that takes a data frame and returns a heatmap plot
makeHeatMap <- function(df, xLabels, yLabels, xTicksNum = 10){
    xLabels <- colnames(df)
    yLabels <- rownames(df)

    if(is.null(xLabels)){
        xLabels <- seq_len(ncol(df))
    }
    if(is.null(yLabels)){
        yLabels <- seq_len(nrow(df))
    }

    df$y <- rownames(df)
    df_long <- reshape2::melt(df, id.vars = "y", variable.name = "x", value.name = "value")
    colnames(df_long) <- c("y", "x", "value")

    ## sort df_long by rownames(fragMatReorderd)
    df_long$x <- factor(df_long$x, levels = xLabels)
    df_long$y <- factor(df_long$y, levels = rev(yLabels))
    ## show 10 time points on x-axis at most
    if (length(xLabels) > xTicksNum){
        step <- ceiling(length(xLabels) / xTicksNum)
        breaksIdx <- seq(1, length(xLabels), by = step)
        breaks <- xLabels[breaksIdx]
    } else {
        breaks <- xLabels
    }

    ggplot2::ggplot(df_long) +
        ggplot2::geom_tile(ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value)) +
        ggplot2::scale_x_discrete(labels = breaks, breaks = breaks) +
        ggplot2::theme(plot.title = ggtext::element_markdown(hjust = 0.5)) +
        viridis::scale_fill_viridis(option = "turbo") +
        ggplot2::theme_minimal()
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
visuIEEGData <- function(epoch) {
    if (is(epoch, "matrix")){
        epoch <- Epoch(epoch)
    }

    gaps <- 2

    elecNames <- epoch$electrodes
    data <- epoch$data
    elecNum <- nrow(data)
    timesNum <- ncol(data)

    plotData <- standardizeIEEG(data)

    times <- epoch$times
    if (is.null(times)) {
        xlabel <- "Time Index"
        timeTicks <- seq_len(timesNum)
    } else {
        xlabel <- "Time (s)"
        timeTicks <- times
    }

    plotData <- apply(plotData, 1, function(x) x - mean(x))
    plotData <- as.data.frame(plotData)
    plotData$timeTicks <- timeTicks
    breakplot <- (seq_len(elecNum) - 1) * gaps
    
    elecNamesReversed <- rev(elecNames)

    ## add gaps between electrodes
    for (i in seq_along(elecNamesReversed)) {
        elec <- elecNamesReversed[i]
        plotData[[elec]] <- plotData[[elec]] + (i-1) * gaps
    }


    p <- ggplot2::ggplot(data = plotData)
    for (i in seq_along(elecNamesReversed)) {
        elec <- elecNamesReversed[i]
        p <- p + ggplot2::geom_line(ggplot2::aes(x = .data$timeTicks, y = .data[[elec]]))
    }

    p +
        ggplot2::labs(x = xlabel, y = "Electrode", size = 2) +
        ggplot2::scale_y_continuous(labels = elecNamesReversed, breaks = breakplot)
}



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
    elecNum <- nrow(fragMat)
    windowNum <- ncol(fragMat)

    elecNames <- frag$electrodes
    sozIndex <- checkIndex(sozIndex, elecNames)

    group1 <- sozIndex
    group2 <- setdiff(seq_len(elecNum), sozIndex)

    elecColor <- rep("blue", elecNum)
    elecColor[seq_along(group2)] <- "black"

    startTime <- frag$startTimes
    if (is.null(startTime)) {
        xlabel <- "Time Index"
        stimes <- seq_len(windowNum)
    } else {
        xlabel <- "Time (s)"
        stimes <- startTime
    }

    rownames(fragMat) <- frag$electrodes
    colnames(fragMat) <- stimes

    ## prepare the data.frame for visualization
    allIndex <- c(group1, group2)
    df <- fragMat[allIndex, ]


    ## make fragDf a long format data.frame
    ## three columns: Electrode, Time, Value
    df_long <- as.data.frame(as.table(fragMatReorderd))
    colnames(df_long) <- c("Electrode", "Time", "Value")

    ## sort df_long by rownames(fragMatReorderd)
    df_long$Electrode <- factor(df_long$Electrode, levels = rev(rownames(fragMatReorderd)))
    df_long$Time <- factor(df_long$Time, levels = stimes)

    ## show 10 time points on x-axis at most
    step <- ceiling(windowNum / 10)
    breaksIdx <- seq(1, length(stimes), by = step)
    breaks <- stimes[breaksIdx]
    xLabels <- stimes[breaksIdx]

    

    makeHeatMap(df) +
        ggplot2::labs(x = xlabel, y = "Electrode", size = 2) +
        ggplot2::theme(
            axis.text.y = ggtext::element_markdown(size = 6, colour = elecColor), # Adjust depending on electrodes
        )
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
plotFragQuantile <- function(frag, sozIndex = NULL) {
    if (is.null(sozIndex)) {
        sozIndex <- estimateSOZ(frag)
    }
    windowNum <- ncol(frag)

    stat <- fragStat(frag, sozIndex)
    qmatrix <- as.data.frame(stat$qmatrix)

    startTimes <- frag$startTimes
    if (is.null(startTimes)) {
        xlabel <- "Time Index"
        timeTicks <- seq_len(windowNum)
    } else {
        xlabel <- "Time (s)"
        timeTicks <- startTimes
    }

    colnames(qmatrix) <- timeTicks

    makeHeatMap(qmatrix)+
        ggplot2::labs(x = xlabel, y = "Quantiles", size = 2) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 4), # Adjust depending on electrodes
        )
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
plotFragDistribution <- function(frag, sozIndex = NULL) {
    if (is.null(sozIndex)) {
        sozIndex <- estimateSOZ(frag)
    }
    
    sozIndex <- checkIndex(sozIndex, frag$electrodes)

    fragMat <- frag$frag
    windowNum <- ncol(fragMat)

    SOZMat <- fragMat[sozIndex, , drop = FALSE]
    RefMat <- fragMat[-sozIndex, , drop = FALSE]
    
    meanSOZ <- apply(fragMat, 2, mean, na.rm = TRUE)
    semSOZ <- apply(fragMat, 2, function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))
    
    meanRef <- apply(RefMat, 2, mean, na.rm = TRUE)
    semRef <- apply(RefMat, 2, function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))

    startTimes <- frag$startTimes
    if (is.null(startTimes)) {
        xlabel <- "Time Index"
        timeTicks <- seq_len(windowNum)
    } else {
        xlabel <- "Time (s)"
        timeTicks <- startTimes
    }

    upperSOZ <- meanSOZ + semSOZ
    lowerSOZ <- meanSOZ - semSOZ
    upperRef <- meanRef + semRef
    lowerRef <- meanRef - semRef

    plotData <- data.frame(
        timeTicks = timeTicks,
        meanSOZ = meanSOZ,
        upperSOZ = upperSOZ,
        lowerSOZ = lowerSOZ,
        meanRef = meanRef,
        upperRef = upperRef,
        lowerRef = lowerRef
    )

    colors <- c("SOZ +/- sem" = "red", "SOZc +/- sem" = "black")
    ggplot2::ggplot(plotData, ggplot2::aes(x = .data$timeTicks)) +
        ggplot2::xlab(xlabel) +
        ggplot2::ylab("Fragility") +
        ggplot2::geom_line(ggplot2::aes(y = .data$meanSOZ, color = "SOZ +/- sem")) +
        ggplot2::geom_line(ggplot2::aes(y = .data$upperSOZ), color = "red", linetype = "dotted") +
        ggplot2::geom_line(ggplot2::aes(y = .data$lowerSOZ), color = "red", linetype = "dotted") +
        ggplot2::geom_line(ggplot2::aes(y = .data$meanRef, color = "SOZc +/- sem")) +
        ggplot2::geom_line(ggplot2::aes(y = .data$upperRef), color = "black", linetype = "dotted") +
        ggplot2::geom_line(ggplot2::aes(y = .data$lowerRef), color = "black", linetype = "dotted") +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lowerSOZ, ymax = .data$upperSOZ), fill = "red", alpha = 0.5) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lowerRef, ymax = .data$upperRef), fill = "black", alpha = 0.5) +
        ggplot2::scale_color_manual(name = "Electrode groups", values = c(colors))
}
