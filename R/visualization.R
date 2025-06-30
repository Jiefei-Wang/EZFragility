# A plot function that takes a data frame and returns a heatmap plot
makeHeatMap <- function(df, xTicksNum = 10, maxLabels = Inf){
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

    ## limit the number of labels on y-axis
    yLabelsForDisplay <- rev(yLabels)  # Match the reversed factor levels
    if (length(yLabelsForDisplay) > maxLabels) {
        by_num <- ceiling(length(yLabelsForDisplay)/maxLabels)
        label_idx <- seq(length(yLabelsForDisplay), 1, by=-by_num)
        yLabelsForDisplay[-label_idx] <- ""
    }

    ggplot(df_long) +
        geom_tile(aes(x = .data$x, y = .data$y, fill = .data$value)) +
        scale_x_discrete(labels = breaks, breaks = breaks) +
        scale_y_discrete(labels = yLabelsForDisplay, breaks = rev(yLabels)) +
        theme(plot.title = element_markdown(hjust = 0.5)) +
        scale_fill_viridis(option = "turbo") +
        theme_minimal()
}


#' Visualization of ictal iEEG
#'
#' @inheritParams calcAdjFrag
#' @inheritParams fragStat
#' @param maxLabels Integer. Maximum number of labels to show on y-axis. Default is 50. The actual number of labels may be less than this value if there are too many electrodes.
#' @return A ggplot object
#'
#' @examples
#' data("pt01EcoG")
#' 
#' ## Visualize a subject of electrodes
#' sozIndex <- which(rowData(pt01EcoG)$soz)
#' display <- c(sozIndex, 77:80)
#' 
#' visuIEEGData(epoch = pt01EcoG[display, ])
#' @export
visuIEEGData <- function(epoch, groupIndex = NULL, maxLabels = 50) {
    if (is(epoch, "matrix")){
        epoch <- Epoch(epoch)
    }

    gaps <- 2

    elecNames <- rownames(epoch)
    data <- tblData(epoch)
    elecNum <- nrow(data)
    timesNum <- ncol(data)

    # group electrodes
    groupIndex <- checkIndex(groupIndex, elecNames)
    group1 <- groupIndex
    group2 <- setdiff(seq_len(elecNum), groupIndex)
    
    # reorder the electrodes
    plotData <- data[c(group1, group2), , drop = FALSE]
    elecNames <- c(elecNames[group1], elecNames[group2])
    plotData <- standardizeIEEG(plotData)

    timePoints <- coltimes(epoch)
    if (is.null(timePoints)) {
        xlabel <- "Time Index"
        timeTicks <- seq_len(timesNum)
    } else {
        xlabel <- "Time (s)"
        timeTicks <- timePoints
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

    elecColor <- rep("blue", elecNum)
    elecColor[seq_along(group2)] <- "black"

    p <- ggplot(data = plotData)
    for (i in seq_along(elecNamesReversed)) {
        elec <- elecNamesReversed[i]
        p <- p + geom_line(aes(x = .data$timeTicks, y = .data[[elec]]))
    }

    ## limit the number of labels on y-axis
    if (length(elecNamesReversed) > maxLabels) {
        by_num <- ceiling(length(elecNamesReversed)/maxLabels)
        label_idx <- seq(length(elecNamesReversed), 1, by=-by_num)
        elecNamesReversed[-label_idx] <- ""
    }

    p +
        labs(x = xlabel, y = "Electrode", size = 2) +
        scale_y_continuous(labels = elecNamesReversed, breaks = breakplot) +
        theme(
            axis.text.y = element_markdown(colour = elecColor)
        )
}



#' Visualization functions (raw signal, fragility matrix)
#'
#' @description `plotFragHeatmap`: plot fragility heatmaps with electrodes marked as soz colored
#'
#' @param frag Fragility object from \code{calcAdjFrag}
#' @inheritParams fragStat
#' @inheritParams visuIEEGData
#' 
#' @return A ggplot object
#'
#' @examples
#' 
#' data("pt01EcoG")
#' 
#' ## sozNames is the name of the electrodes we assume are in the SOZ
#' sozNames <- metaData(pt01EcoG)$sozNames
#' 
#' ## precomputed fragility object
#' data("pt01Frag")
#' 
#' ## plot the fragility heatmap
#' plotFragHeatmap(frag = pt01Frag, groupIndex = sozNames)
#' 
#' @rdname plotFragHeatmap
#' @export
plotFragHeatmap <- function(
    frag,
    groupIndex = NULL,
    maxLabels = 50) {
    fragMat <- frag$frag
    elecNum <- nrow(fragMat)
    windowNum <- ncol(fragMat)

    elecNames <- frag$electrodes
    groupIndex <- checkIndex(groupIndex, elecNames)

    group1 <- groupIndex
    group2 <- setdiff(seq_len(elecNum), groupIndex)

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

    makeHeatMap(df, maxLabels = maxLabels) +
        labs(x = xlabel, y = "Electrode", size = 2) +
        theme(
            axis.text.y = element_markdown(size = 6, colour = elecColor), # Adjust depending on electrodes
        )
}


#' @description `plotFragQuantile`: Plot Fragility time quantiles for two electrodes groups
#' 
#' @rdname plotFragHeatmap
#' @examples
#' ## plot the fragility quantiles
#' plotFragQuantile(frag = pt01Frag, groupIndex = sozNames)
#' 
#' @export
plotFragQuantile <- function(frag, groupIndex = NULL, groupName = "SOZ") {
    if (is.null(groupIndex)) {
        groupIndex <- estimateSOZ(frag)
    }
    groupIndex <- checkIndex(groupIndex, frag$electrodes)
    windowNum <- ncol(frag)

    stat <- fragStat(
        frag, 
        groupIndex = groupIndex, 
        groupName = groupName
    )
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
        labs(x = xlabel, y = "Quantiles", size = 8) +
        theme(
            axis.text.y = element_text(size = 8), # Adjust depending on electrodes
        )
}


#' @description `plotFragQuantile`: Plot Fragility time distribution for two electrodes groups
#' @param bandType Character. The type of band to use, either "SEM" or "SD". Default is "SEM".
#' @rdname plotFragHeatmap
#' @examples
#' ## plot the fragility distribution
#' plotFragDistribution(frag = pt01Frag, groupIndex = sozNames)
#' 
#' @export
plotFragDistribution <- function(frag, groupIndex = NULL, groupName="SOZ", bandType = c("SEM", "SD")) {
    bandType <- match.arg(bandType)
    if (is.null(groupIndex)) {
        groupIndex <- estimateSOZ(frag)
    }
    

    windowNum <- ncol(frag$frag)
    stat <- fragStat(
        frag, 
        groupIndex = groupIndex
    )

    groupMean <- stat$groupMean
    refMean <- stat$refMean
    if (bandType == "SEM") {
        groupWidth <- stat$groupSEM
        refWidth <- stat$refSEM
    } else if (bandType == "SD") {
        groupWidth <- stat$groupSD
        refWidth <- stat$refSD
    }
    groupUpperBound <- groupMean + groupWidth
    groupLowerBound <- groupMean - groupWidth
    refUpperBound <- refMean + refWidth
    refLowerBound <- refMean - refWidth

    startTimes <- frag$startTimes
    if (is.null(startTimes)) {
        xlabel <- "Time Index"
        timeTicks <- seq_len(windowNum)
    } else {
        xlabel <- "Time (s)"
        timeTicks <- startTimes
    }


    plotData <- data.frame(
        timeTicks = timeTicks,
        groupMean = groupMean,
        groupUpperBound = groupUpperBound,
        groupLowerBound = groupLowerBound,
        refMean = refMean,
        refUpperBound = refUpperBound,
        refLowerBound = refLowerBound
    )
    
    groupColor <- glue("Group +/- {bandType}")
    refColor <- glue("Ref +/- {bandType}")
    colors <- setNames(c("red", "black"), c(groupColor, refColor))

    ggplot(plotData, aes(x = .data$timeTicks)) +
        xlab(xlabel) +
        ylab("Fragility") +
        geom_line(
            aes(y = .data$groupMean, color = groupColor)
        ) +
        geom_line(aes(y = .data$groupUpperBound), color = "red", linetype = "dotted") +
        geom_line(aes(y = .data$groupLowerBound), color = "red", linetype = "dotted") +
        geom_line(aes(y = .data$refMean, color = refColor)) +
        geom_line(aes(y = .data$refUpperBound), color = "black", linetype = "dotted") +
        geom_line(aes(y = .data$refLowerBound), color = "black", linetype = "dotted") +
        geom_ribbon(aes(ymin = .data$groupLowerBound, ymax = .data$groupUpperBound), fill = "red", alpha = 0.5) +
        geom_ribbon(aes(ymin = .data$refLowerBound, ymax = .data$refUpperBound), fill = "black", alpha = 0.5) +
        scale_color_manual(name = "Electrode groups", values = c(colors))
}
