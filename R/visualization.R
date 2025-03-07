#' Visualization functions (raw signal, fragility matrix)
#' 
#' plot fragility heatmaps with electrodes marked as soz colored
#' 
#' @inheritParams fragStat
#' @param sozID Integer or string. Vector soz electrodes 
#' @param timeRange Numeric Vector of length 2. The time window to display at the x-axis
#' @param title String. Figure title
#' @param display Integer or string. Vector electrodes to display
#'
#' @return Heatmap plot of the fragility matrix with soz electrodes in blue in the bottom
#'
#' @examples
#'# use integer index for display and soz electrodes
#'data("pt01Epochm1sp2s")
#'sozIndex<-attr(pt01Epochm1sp2s,"sozIndex")
#'data("pt01Fragm1sp2s")
#'timeRange <- c(-1,2)
#'display <- c(sozIndex,77:80)
#'fragplot<-heatmapFrag(frag=pt01Fragm1sp2s,sozID=sozIndex,
#'timeRange = timeRange,title="PT01 seizure 1",display=display)
#'fragplot
#'
#'
#'# use electrodes name for display and soz electrodes
#'data("pt01Epochm1sp2s")
#'sozNames<-attr(pt01Epochm1sp2s,"sozNames")
#'data("pt01Fragm1sp2s")
#'timeRange <- c(-1,2)
#'display <- c(sozNames,"MLT1","MLT2","MLT3","MLT4")
#'fragplot<-heatmapFrag(frag=pt01Fragm1sp2s,sozID=sozNames,
#'                     timeRange = timeRange,title="PT01 seizure 1",display=display)
#'fragplot
#'
#' # save plot to file with ggplot2
#'data("pt01Epochm1sp2s")
#'data("pt01Fragm1sp2s")
#'sozIndex<-attr(pt01Epochm1sp2s,"sozIndex")
#'timeRange <- c(-10,10)
#'display <- c(sozIndex,77:80)
#'pathplot <- "~"
#'title <- "PT01sz1"
#'resfile <- paste(pathplot,'/FragilityHeatMap',title,'.png',sep="")
#'fragplot<-heatmapFrag(frag=pt01Fragm1sp2s,sozID=sozIndex,timeRange=timeRange,
#'                      title=title,display=display)
#'fragplot
#'ggplot2::ggsave(resfile)
#' 
#' @export
heatmapFrag<-function(frag,sozID,timeRange = NULL,title="Patient name seizure number",display=NULL){
  titlepng<-title
  
  
  if(is.null(display)){
    display<-1:nrow(frag$frag)
  }
  
  if (is(frag, "Fragility")) {
    frag <- frag$frag
  }

  
  elecName<-rownames(frag)
  elecInd=c(1:nrow(frag))
  
  
  if(typeof(display)=="integer"){  
    
    displayTot<-1:nrow(frag)
    diffDisplayTot<-setdiff(display,displayTot)
    
    if(length(diffDisplayTot)!=0){
      listDisplayMissing<-paste(as.character(diffDisplayTot),collapse=" ")
      message<-paste("Number(s) ",listDisplayMissing,"are out of electrode number limit")
      warning(message)
      display<-display[!display%in%diffDisplayTot]
      displayCor<-paste(as.character(display),collapse=" ")
      message<-paste("Keeping indices.",displayCor)
      warning(message)
      
    }  
    displayid<-display
    
  }else{
    
    diffDisplayTot<-setdiff(display,elecName)
    
    if(length(diffDisplayTot)!=0){
      listDisplayMissing<-paste(diffDisplayTot,collapse=" ")
      message<-paste(" Name(s) ",listDisplayMissing,"are out of name list")
      warning(message)
      display<-display[!display%in%diffDisplayTot]
      displayCor<-paste(display,collapse=" ")
      message<-paste("Keeping names.",displayCor)
      warning(message)
      
    }  
    
    displayid<-which(elecName%in%display)
    
  }
  
  fragDisplay<-frag[displayid,]
  nElec <- nrow(fragDisplay)
  elecTot<-c(1:nElec)
  
  if(typeof(sozID)=="integer"){  
    
    diffelecInd<-setdiff(sozID,elecInd)
    
    if(length(diffelecInd)!=0){
      listElecMissing<-paste(as.character(diffelecInd),collapse=" ")
      message<-paste("Number(s) ",listElecMissing,"are out of electrode number limit")
      warning(message)
      sozID<-sozID[!sozID%in%diffelecInd]
      listSozCor<-paste(as.character(sozID),collapse=" ")
      message<-paste("Keeping indices.",listSozCor)
      warning(message)
    }
    sozIDid<-sozID
    
  }else{

    diffsoztot<-setdiff(sozID,elecName)
    
    if(length(diffsoztot)!=0){
      listSozMissing<-paste(diffsoztot,collapse=" ")
      message<-paste("Name(s) ",listSozMissing,"are out of name list")
      warning(message)
      sozID<-sozID[!sozID%in%diffsoztot]
      sozcor<-paste(sozID,collapse=" ")
      message<-paste("Keeping names.",sozcor)
      warning(message)
      
    }  
    
    sozIDid<-which(elecName%in%sozID)
  }

  sozIDd<-which(displayid%in%sozIDid)
  sozIDCd<-which(!displayid%in%sozIDid)
  sozIDSozC<-c(sozIDd,sozIDCd)

  elecNum <- rownames(fragDisplay)
  nw<- ncol(fragDisplay)
  colorelec<-elecNum
  nsoz<-length(sozID)
  colorelec[1:nElec]<-"blue"
  nb<-nElec-nsoz
  colorelec[1:nb]<-"black"

  sozIDSozC<-rev(sozIDSozC)
  fragord<-fragDisplay[sozIDSozC,]
  fragdf<-data.frame(fragord)
  
  if(is.null(timeRange)){
    xlabel<-"Time Index"
    stimes<-seq_len(nw)
  }
  else{
    xlabel<-"Time (s)"
    stimes<-seq(timeRange[1],timeRange[2],length.out=nw)
  }

  colnames(fragdf)<-stimes
  rownames(fragdf)<-elecNum
  
  elecNum<-rev(elecNum)
  
  fragmapData <- expand.grid(Time = stimes, Electrode = elecNum)
  fragmapData$Value <- c(t(fragord))

  
  p<-ggplot2::ggplot(fragmapData, ggplot2::aes(x = .data$Time, y = .data$Electrode, fill = .data$Value)) +
    ggplot2::geom_tile() +
    ggplot2::ggtitle(as.character(titlepng)) +
    ggplot2::theme(plot.title=ggtext::element_markdown(hjust=0.5)) +
    ggplot2::labs(x = xlabel, y = "Electrode",size=2) +
    viridis::scale_fill_viridis(option = "turbo") +
    ggplot2::geom_vline(xintercept =0, 
                  color = "black", linetype = "dashed", linewidth = 1) +
   ggplot2::theme_minimal() +
   ggplot2::theme(
      axis.text.y = ggtext::element_markdown(size=6,colour=colorelec),     # Adjust depending on electrodes
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
#'data("pt01Epochm1sp2s")
#'sozIndex <- attr(pt01Epochm1sp2s,"sozIndex")
#'display <- c(sozIndex,77:80)
#'timeRange <- c(-1,2)
#'iEEGplot<-visuIEEGData(ieegts=pt01Epochm1sp2s,timeRange=timeRange,display=display)
#'iEEGplot
#' @export
visuIEEGData<-function(ieegts, timeRange=NULL, title = "Patient name seizure number", display=NULL){
 
  titlepng<- title
  if(is.null(display)){
    display<-1:ncol(ieegts)
  }

  
  elecName<-colnames(ieegts)
  if(typeof(display)=="integer"){  

    displayTot<-1:ncol(ieegts)
    diffDisplayTot<-setdiff(display,displayTot)
    
    if(length(diffDisplayTot)!=0){
      listDisplayMissing<-paste(as.character(diffDisplayTot),collapse=" ")
      message<-paste("Display electrodes indices. Numbers ",listDisplayMissing,"are out of electrode number limit")
      warning(message)
      display<-display[!display%in%diffDisplayTot]
      displayCor<-paste(as.character(display),collapse=" ")
      message<-paste("Keeping indices.",displayCor)
      warning(message)
      
    }  
  
     displayid<-display
     
  }else{

    diffDisplayTot<-setdiff(display,elecName)
    
    if(length(diffDisplayTot)!=0){
      listDisplayMissing<-paste(diffDisplayTot,collapse=" ")
      message<-paste("Display electrodes names. Names ",listDisplayMissing,"are out of name list")
      warning(message)
      display<-display[!display%in%diffDisplayTot]
      displayCor<-paste(display,collapse=" ")
      message<-paste("Keeping names.",displayCor)
      warning(message)
      
    }  
        
    displayid<-which(elecName%in%display)
  }
  
  scaling <- 10^floor(log10(max(ieegts)))
  plotData<-ieegts[,displayid]/scaling
  gaps<-2
  displayNames<-colnames(ieegts)[displayid]
  nElec<-length(displayid)
  nt<-nrow(plotData)
  if(is.null(timeRange)){
    xlabel<-"Time Index"
    stimes<-seq_len(nt)
  }else{
    xlabel<-"Time (s)"
    stimes<-seq(timeRange[1],timeRange[2],length.out=nt)
  }

  
  for(i in 1:ncol(plotData)){
    plotData[, i] <- (plotData[, i]- mean(plotData[, i]))+
      (ncol(plotData)-i)*gaps
  }
  plotData<-data.frame(plotData)
  breakplot<-(c(1:nElec)-1)*gaps
  
  p<-ggplot2::ggplot(data=plotData,ggplot2::aes(x=stimes,y=plotData))+
    ggplot2::ggtitle(titlepng)+
    ggplot2::labs(x = xlabel, y = "Electrode",size=2)+ 
    ggplot2::geom_vline(xintercept =0, 
                        color = "black", linetype = "dashed", size = 1)
  
  for(i in 1:nElec){
    p<-p+ggplot2::geom_line(ggplot2::aes_string(y=names(plotData)[i]))
  }
  displayNames<-rev(displayNames)
  p<-p+ggplot2::scale_y_continuous(labels=displayNames,breaks=breakplot)
  
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
#'timeRange <- c(-1,2)
#'data("pt01Epochm1sp2s")
#'sozIndex<-attr(pt01Epochm1sp2s,"sozIndex")
#'data("pt01Fragm1sp2s")
#'# compute fragility statistics evolution with time (mean and standard deviation) for soz and
#'# non soz groups
#'pt01fragstat <- fragStat(frag=pt01Fragm1sp2s, sozID=sozIndex)
#'plotFragQuantile(FragStatObj=pt01fragstat, timeRange=timeRange)
plotFragQuantile<-function(FragStatObj, timeRange = NULL,title="Fragility Quantiles over time"){
  if(is(FragStatObj, "FragStat")){
    qmatrix <- FragStatObj$qmatrix
  }
  
  nw <- ncol(qmatrix)
  if(is.null(timeRange)){
    xlabel<-"Time Index"
    stimes<-seq_len(nw)
  }
  else{
    xlabel<-"Time (s)"
    stimes<-seq(timeRange[1],timeRange[2],length.out=nw)
  }


  quantilesName<-rownames(qmatrix)
  quantilePlot<- expand.grid(Time = stimes, Stats=quantilesName)
  quantilePlot$Value <- c(t(qmatrix))
  
  titlepng <- title
  
  p<-ggplot2::ggplot(quantilePlot, ggplot2::aes(x = .data$Time, y = .data$Stats, fill = .data$Value)) +
    ggplot2::geom_tile() +
    ggplot2::ggtitle(titlepng)+
    ggplot2::labs(x = xlabel, y = "Quantiles",size=2) +
    viridis::scale_fill_viridis(option = "turbo") +  #
    
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size=4),     # Adjust depending on electrodes
    )

  if(!is.null(timeRange)){
    p <- p + ggplot2::geom_vline(xintercept =0, 
                                 color = "black", linetype = "dashed", linewidth = 1)}

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
#'data("pt01Epochm1sp2s")
#'sozIndex<-attr(pt01Epochm1sp2s,"sozIndex")
#'# Load the precomputed fragility object
#'timeRange <- c(-10,10)
#'data("pt01Fragm1sp2s")
#'# compute fragility statistics evolution with time (mean and standard deviation) for soz and
#'# non soz groups
#'pt01fragstat <- fragStat(frag=pt01Fragm1sp2s, sozID=sozIndex)
#'# plot the statistical results
#'pfragstat<-plotFragDistribution(stat=pt01fragstat,timeRange=timeRange)
#'pfragstat
plotFragDistribution<-function(
  stat = NULL, 
  timeRange = NULL,
  title='Average Fragility over time'){
  stopifnot(is(stat, "FragStat"))
  cmeansoz <- stat$cmeansoz
  cmeansozc <- stat$cmeansozc
  csdsoz <- stat$csdsoz
  csdsozc <- stat$csdsozc
  
  nw <- length(cmeansoz)
  if(is.null(timeRange)){
    xlabel<-"Time Index"
    stimes<-seq_len(nw)
  }
  else{
    xlabel<-"Time (s)"
    stimes<-seq(timeRange[1],timeRange[2],length.out=nw)
  }

  
  sozsdp <- cmeansoz+csdsoz
  sozsdm <- cmeansoz-csdsoz
  sozcsdp <- cmeansozc+csdsozc
  sozcsdm <- cmeansozc-csdsozc
  
  plotmeanstd<-as.data.frame(stimes)
  colnames(plotmeanstd)<-"times"
  plotmeanstd$meansoz<-cmeansoz
  plotmeanstd$sozsdp<-sozsdp
  plotmeanstd$sozsdm<-sozsdm
  plotmeanstd$meansozc<-cmeansozc
  plotmeanstd$sozcsdp<-sozcsdp
  plotmeanstd$sozcsdm<-sozcsdm
  
  titlepng <- title
  colors<-c("SOZ +/- sem" = "red", "SOZc +/- sem" = "black")
  ggplot2::theme_grey(base_size = 22)
  p<-ggplot2::ggplot(plotmeanstd, ggplot2::aes(x=.data$times, y=.data$meansoz))+ 
   ggplot2::xlab(xlabel)+
   ggplot2::ylab('Fragility')+
   ggplot2::ggtitle(titlepng)+
   ggplot2::geom_line(ggplot2::aes(y = .data$meansoz,color="SOZ +/- sem"))+  
   ggplot2::geom_line(ggplot2::aes(y = .data$sozsdp),color='red',linetype="dotted")+  
   ggplot2::geom_line(ggplot2::aes(y = .data$sozsdm),color='red',linetype="dotted")+ 
   ggplot2::geom_line(ggplot2::aes(y = .data$meansozc,color="SOZc +/- sem"))+
   ggplot2::geom_line(ggplot2::aes(y = .data$sozcsdp),color='black',linetype="dotted")+  
   ggplot2::geom_line(ggplot2::aes(y = .data$sozcsdm),color='black',linetype="dotted")+
   ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$sozsdm,ymax=.data$sozsdp), fill="red",alpha=0.5)+
   ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$sozcsdm,ymax=.data$sozcsdp), fill="black",alpha=0.5)+  
   ggplot2::scale_color_manual(name="Electrode groups",values = c(colors)) 
    
  ## add vertical line at time 0 if timeRange is specified
  if(!is.null(timeRange)){
    p <- p + ggplot2::geom_vline(xintercept =0, 
                  color = "black", linetype = "dashed", linewidth = 1)}
  return(p)
}
