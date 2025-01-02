#' Visualization functions (raw signal, fragility matrix)
#' 
#' plot fragility heatmaps with electrodes marked as soz colored
#'
#' @param frag Numeric. Fragility matrix results
#' @param elecsoz Integer. Vector soz electrodes (for good electrodes)
#' @param time_window Numeric Vector. Fragility heatmap time window around seizure onset (s)
#' @param title String. Figure title
#' @param display Integer. Electrodes to display
#'
#' @return Heatmap plot of the fragility matrix with soz electrodes in blue in the bottom
#'
#' @examples
#' data("fragm3sp5s")
#' data("elecsoz")
#' time_window=c(-3:5)
#' display=c(elecsoz,77:80)
#' heatmap_frag(frag=fragm3sp5s,elecsoz=elecsoz,time_window=c(-3,5),display=display)
#' @export
heatmap_frag<-function(frag,elecsoz,time_window,option=NULL,title="Fragility heatmap",display=NULL){
  titlepng<-title
  if(is.null(display)){
    display<-1:nrow(frag)
  }
  fragdisplay<-frag[display,]
  n_elec <- nrow(fragdisplay)
  electot<-c(1:n_elec)
  
  elecsozd=which(display%in%elecsoz)
  elecsozcd=which(!display%in%elecsoz)
  elecsozsozc=c(elecsozd,elecsozcd)

  elecnum <- rownames(fragdisplay)
  nw<- ncol(fragdisplay)
  colorelec<-elecnum
  nsoz=length(elecsoz)
  colorelec[1:n_elec]="black"
  colorelec[1:nsoz]="blue"

  fragord<-fragdisplay[elecsozsozc,]
  fragdf<-data.frame(fragord)
  stimes=c(1:nw)*(time_window[2]-time_window[1])/nw+time_window[1]
  colnames(fragdf)<-stimes
  rownames(fragdf)<-elecnum
  
  fragmap_data <- expand.grid(Time = stimes, Electrode = elecnum)
  fragmap_data$Value <- c(t(fragord))
  
  ggplot2::ggplot(fragmap_data, ggplot2::aes(x = Time, y = Electrode, fill = Value)) +
    ggplot2::geom_tile() +
    ggplot2::ggtitle(titlepng)+
    ggplot2::labs(x = "Time (s)", y = "Electrode",size=2) +
    viridis::scale_fill_viridis(option = "turbo") +  #
    
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size=4,colour=colorelec),     # Adjust depending on electrodes
    )
}

#' Visualization of ictal iEEG 
#'
#' @param ieegts Numeric. A matrix of iEEG time series x(t),
#' with time points as rows and electrodes names as columns
#' @param scaling Numeric. Scaling factor
#' @param displayChannels Integer. Channels evctor to display
#'
#' @return plot raw signal
#'
#' @examples
#' data("PT01Epochm30sp30s")
#' data("ElectrodesDataPT01")
#' displayChannels=which(ElectrodesDataPT01$insoz==TRUE)
#' visuiEEGdata(ieegts=PT01Epochm30sp30s,1000000, displayChannels = displayChannels)
#' @export
visuiEEGdata<-function(ieegts, scaling, displayChannels){
  plotData<-ieegts[,displayChannels]/scaling
  gaps<-2
  displayNames=colnames(ieegts)[displayChannels]
  
  for(i in seq_along(plotData)){
    plotData[, i] <- (plotData[, i]- mean(plotData[, i]))+
      (ncol(plotData)-i)*gaps
  }
  
  plot(plotData[, 1],type="l" ,cex=0.1,
       ylim = range(plotData), yaxt = "n")
  for(i in 2:ncol(plotData)){
    lines(plotData[, i])
  }
  axis(2, at = rev(seq_along(displayChannels) - 1)*gaps,
       labels = displayNames,las=1)
  
  
}

#' Plot Fragility time quantiles for two electrodes group marked as soz non marked as soz
#'
#' @param qmatsozsozc Numeric. quantile matrix for the two groups
#' @param time_window_ictal Numeric. time window of processed ictal iEEG
#' @param title String. Figure title
#'
#' @return Quantile plot
#' @export
#'
#' @examples
#' data("fragstat")
#' time_window_ictal=c(-3,5)
#' plot_frag_quantile( qmatsozsozc=fragstat[[1]], time_window_ictal=time_window_ictal,title=title)
plot_frag_quantile<-function( qmatsozsozc, time_window_ictal,title=title){
 
  nw=ncol(qmatsozsozc)
  stimes=c(1:nw)*(time_window_ictal[2]-time_window_ictal[1])/nw+time_window_ictal[1] 
  quantilesname<-c(paste0("SOZ(", seq(10,100,by=10),")"),paste0("SOZc(", seq(10,100,by=10),")"))
  quantileplot<- expand.grid(Time = stimes, Stats=quantilesname)
  quantileplot$Value <- c(t(qmatsozsozc))
  
  titlepng=title
  
  ggplot2::ggplot(quantileplot, ggplot2::aes(x = Time, y = Stats, fill = Value)) +
    ggplot2::geom_tile() +
    ggplot2::ggtitle(titlepng)+
    ggplot2::labs(x = "Time (s)", y = "Quantiles",size=2) +
    viridis::scale_fill_viridis(option = "turbo") +  #
    
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size=4),     # Adjust depending on electrodes
    )
  
 
}


#'  Plot Fragility time distribution for two electrodes group marked and non-marked as soz
#'
#' @param statcurves. List of four numeric vectors for mean and standard deviation statistics of the two groups
#' @param time_window_ictal time window of processed ictal iEEG
#'
#' @return plot fragility distribution
#' @export
#'
#' @examples
#' time_window_ictal=c(-3,5)
#' fragstatcurve=fragstat[-1]
#' plot_frag_distribution(statcurves=fragstatcurve,time_window_ictal=time_window_ictal,title="PT01 seizure 1")

plot_frag_distribution<-function( statcurves,time_window_ictal,title=title){
  
  cmeansoz=statcurves[[1]] # mean soz group in function of time window
  cmeansozc=statcurves[[2]] # mean non soz group in function of time window
  csdsoz=statcurves[[3]] # standard deviation soz group in function of time window
  csdsozc=statcurves[[4]]  # standard deviation non soz group in function of time window
  
  nw=length(cmeansoz)
  stimes=c(1:nw)*(time_window_ictal[2]-time_window_ictal[1])/nw+time_window_ictal[1] 
  
  sozsdp=cmeansoz+csdsoz
  sozsdm=cmeansoz-csdsoz
  sozcsdp=cmeansozc+csdsozc
  sozcsdm=cmeansozc-csdsozc
  
  plotmeanstd<-as.data.frame(stimes)
  colnames(plotmeanstd)<-"times"
  plotmeanstd$meansoz<-cmeansoz
  plotmeanstd$sozsdp<-sozsdp
  plotmeanstd$sozsdm<-sozsdm
  plotmeanstd$meansozc<-cmeansozc
  plotmeanstd$sozcsdp<-sozcsdp
  plotmeanstd$sozcsdm<-sozcsdm
  
  titlepng=title
  
  ggplot2::ggplot(plotmeanstd, ggplot2::aes(x=times, y=cmeansoz)) + 
    ggplot2::xlab('Time around seizure in s')+
    ggplot2::ylab('Fragility')+
    ggplot2::ggtitle(titlepng)+
    ggplot2::geom_line(ggplot2::aes(y = meansoz),color='red') + 
    ggplot2::geom_line(ggplot2::aes(y = sozsdp),color='red',linetype="dotted") + 
    ggplot2::geom_line(ggplot2::aes(y = sozsdm),color='red',linetype="dotted") + 
    ggplot2::geom_line(ggplot2::aes(y = meansozc),color='black')+ 
    ggplot2::geom_line(ggplot2::aes(y = sozcsdp),color='black',linetype="dotted") + 
    ggplot2::geom_line(ggplot2::aes(y = sozcsdm),color='black',linetype="dotted")+ 
    ggplot2::geom_ribbon(ggplot2::aes(ymin=sozsdm,ymax=sozsdp), fill="red",alpha=0.5)+
    ggplot2::geom_ribbon(ggplot2::aes(ymin=sozcsdm,ymax=sozcsdp), fill="black",alpha=0.5)  
  

}
