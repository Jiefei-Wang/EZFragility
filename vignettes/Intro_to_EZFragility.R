## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(EZFragility)

## ----ictal data---------------------------------------------------------------
data("pt01Epochm1sp2s")
# seizure onset (soz) electrodes indices
sozIndex <- attr(pt01Epochm1sp2s,"sozIndex")
# seizure onset electrodes names
sozNames <- attr(pt01Epochm1sp2s,"sozZames")

## Visualize ictal data
display <- c(sozIndex,77:80)
timeRange <- c(-1,2)
iEEGplot<-visuIEEGData(ieegts=pt01Epochm1sp2s,timeRange=timeRange,display=display)
iEEGplot

## ----compute fragility, eval=FALSE--------------------------------------------
# window <- 250
# step <- 125
# lambda <- NULL
# nSearch <- 10
# pt01Fragm1sp2s <- calcAdjFrag(ieegts = pt01Epochm1sp2s, window = window, step = step, lambda = lambda,nSearch=nSearch)
# # Fragility matrix result
# pt01Fragm1sp2s$frag[1:5, 1:5]

## -----------------------------------------------------------------------------
data("pt01Epochm1sp2s")
sozIndex<-attr(pt01Epochm1sp2s,"sozIndex")
data("pt01Fragm1sp2s")
timeRange <- c(-1,2)
display <- c(sozIndex,77:80)
fragplot<-heatmapFrag(frag=pt01Fragm1sp2s,sozID=sozIndex,timeRange = timeRange,title="PT01 seizure 1",display=display)
fragplot

## ----compute_frag_stat, eval=TRUE---------------------------------------------
fragstat <- fragStat(frag=pt01Fragm1sp2s, sozID=sozIndex)
fragstat

## ----visualize_frag_auantile,out.width="100%"---------------------------------
# plot the statistical results
pfragstat<-plotFragDistribution(stat =fragstat,timeRange=timeRange)
pfragstat

## ----visualize_frag_quantiles,out.width="100%"--------------------------------
qfragplot<-plotFragQuantile(FragStatObj=fragstat, timeRange=timeRange)
qfragplot

