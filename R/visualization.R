# A plot function that takes a data frame and returns a heatmap plot
makeHeatMap <- function(df, xTicksNum = 10){
    xLabels <- colnames(df)
    yLabels <- rownames(df)

    if(is.null(xLabels)){
        xLabels <- seq_len(ncol(df))
    }
    if(is.null(yLabels)){
        yLabels <- seq_len(nrow(df))
    }

    df$y <- yLabels
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
#' @param epoch Matrix or Epoch object. iEEG data matrix or Epoch object. If matrix, the row names are the electrode names and the column names are the time points
#' @return A ggplot object
#'
#' @examples
#' data("pt01EcoG")
#' 
#' ## Visualize a subject of electrodes
#' sozIndex <- attr(pt01EcoG, "sozIndex")
#' display <- c(sozIndex, 77:80)
#' 
#' epoch <- Epoch(pt01EcoG)
#' visuIEEGData(epoch = epoch[display, ])
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
#' @description `plotFragHeatmap`: plot fragility heatmaps with electrodes marked as soz colored
#'
#' @param frag Fragility object from \code{calcAdjFrag}
#' @param sozIndex Integer or string. A group of electrodes to mark as in the Seizure Onset Zone (SOZ)
#' 
#' @return A ggplot object
#'
#' @examples
#' 
#' data("pt01EcoG")
#' 
#' ## sozIndex is the index of the electrodes we assume are in the SOZ
#' sozIndex <- attr(pt01EcoG, "sozIndex")
#' 
#' ## index of the electrodes to display
#' display <- c(sozIndex, 77:80)
#' 
#' ## precomputed fragility object
#' data("pt01Frag")
#' 
#' ## plot the fragility heatmap
#' plotFragHeatmap(frag = pt01Frag[display, ], sozIndex = sozIndex)
#' 
#' @rdname plotFragHeatmap
#' @export
plotFragHeatmap <- function(
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
    df <- as.data.frame(fragMat[allIndex, ])



    makeHeatMap(df) +
        ggplot2::labs(x = xlabel, y = "Electrode", size = 2) +
        ggplot2::theme(
            axis.text.y = ggtext::element_markdown(size = 6, colour = elecColor), # Adjust depending on electrodes
        )
}


#' @description `plotFragQuantile`: Plot Fragility time quantiles for two electrodes group marked as SOZ and reference
#' 
#' @rdname plotFragHeatmap
#' @examples
#' ## plot the fragility quantiles
#' plotFragQuantile(frag = pt01Frag[display, ], sozIndex = sozIndex)
#' 
#' @export
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


#' @description `plotFragQuantile`: Plot Fragility time distribution for two electrodes group marked as SOZ and reference
#' 
#' @rdname plotFragHeatmap
#' @examples
#' ## plot the fragility distribution
#' plotFragDistribution(frag = pt01Frag[display, ], sozIndex = sozIndex)
#' 
#' @export
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
