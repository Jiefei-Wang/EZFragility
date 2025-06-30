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

    ggplot(df_long) +
        geom_tile(aes(x = .data$x, y = .data$y, fill = .data$value)) +
        scale_x_discrete(labels = breaks, breaks = breaks) +
        theme(plot.title = element_markdown(hjust = 0.5)) +
        scale_fill_viridis(option = "turbo") +
        theme_minimal()
}


#' Visualization of ictal iEEG
#'
#' @inheritParams calcAdjFrag
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
visuIEEGData <- function(epoch) {
    if (is(epoch, "matrix")){
        epoch <- Epoch(epoch)
    }

    gaps <- 2

    elecNames <- rownames(epoch)
    data <- tblData(epoch)
    elecNum <- nrow(data)
    timesNum <- ncol(data)

    plotData <- standardizeIEEG(data)

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


    p <- ggplot(data = plotData)
    for (i in seq_along(elecNamesReversed)) {
        elec <- elecNamesReversed[i]
        p <- p + geom_line(aes(x = .data$timeTicks, y = .data[[elec]]))
    }

    p +
        labs(x = xlabel, y = "Electrode", size = 2) +
        scale_y_continuous(labels = elecNamesReversed, breaks = breakplot)
}



#' Visualization functions (raw signal, fragility matrix)
#'
#' @description `plotFragHeatmap`: plot fragility heatmaps with electrodes marked as soz colored
#'
#' @param frag Fragility object from \code{calcAdjFrag}
#' @param groupIndex Integer or string. A group of electrodes to mark 
#' @param groupName Character. Name of the group of electrodes, default is "SOZ"
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
    groupIndex = NULL) {
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



    makeHeatMap(df) +
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
#' 
#' @rdname plotFragHeatmap
#' @examples
#' ## plot the fragility distribution
#' plotFragDistribution(frag = pt01Frag, groupIndex = sozNames)
#' 
#' @export
plotFragDistribution <- function(frag, groupIndex = NULL, groupName="SOZ") {
    if (is.null(groupIndex)) {
        groupIndex <- estimateSOZ(frag)
    }
    
    groupIndex <- checkIndex(groupIndex, frag$electrodes)

    fragMat <- frag$frag
    windowNum <- ncol(fragMat)

    groupMat <- fragMat[groupIndex, , drop = FALSE]
    RefMat <- fragMat[-groupIndex, , drop = FALSE]
    
    groupMean <- apply(groupMat, 2, mean, na.rm = TRUE)
    groupSEM <- apply(groupMat, 2, function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))

    refMean <- apply(RefMat, 2, mean, na.rm = TRUE)
    refSEM <- apply(RefMat, 2, function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))

    startTimes <- frag$startTimes
    if (is.null(startTimes)) {
        xlabel <- "Time Index"
        timeTicks <- seq_len(windowNum)
    } else {
        xlabel <- "Time (s)"
        timeTicks <- startTimes
    }

    groupUpperBound <- groupMean + groupSEM
    groupLowerBound <- groupMean - groupSEM
    refUpperBound <- refMean + refSEM
    refLowerBound <- refMean - refSEM

    plotData <- data.frame(
        timeTicks = timeTicks,
        groupMean = groupMean,
        groupUpperBound = groupUpperBound,
        groupLowerBound = groupLowerBound,
        refMean = refMean,
        refUpperBound = refUpperBound,
        refLowerBound = refLowerBound
    )
    
    color1 <- "Group +/- sem"
    color2 <- "Groupc +/- sem"
    colors <- c("Group +/- sem" = "red", "Groupc +/- sem" = "black")
    ggplot(plotData, aes(x = .data$timeTicks)) +
        xlab(xlabel) +
        ylab("Fragility") +
        geom_line(
            aes(y = .data$groupMean, color = "Group +/- sem")
        ) +
        geom_line(aes(y = .data$groupUpperBound), color = "red", linetype = "dotted") +
        geom_line(aes(y = .data$groupLowerBound), color = "red", linetype = "dotted") +
        geom_line(aes(y = .data$refMean, color = "Groupc +/- sem")) +
        geom_line(aes(y = .data$refUpperBound), color = "black", linetype = "dotted") +
        geom_line(aes(y = .data$refLowerBound), color = "black", linetype = "dotted") +
        geom_ribbon(aes(ymin = .data$groupLowerBound, ymax = .data$groupUpperBound), fill = "red", alpha = 0.5) +
        geom_ribbon(aes(ymin = .data$refLowerBound, ymax = .data$refUpperBound), fill = "black", alpha = 0.5) +
        scale_color_manual(name = "Electrode groups", values = c(colors))
}
