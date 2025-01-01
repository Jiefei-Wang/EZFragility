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
#' @export
#'
#' @examples
#' data("fragm3sp5s")
#' data("elecsoz")
#' time_window=c(-3:5)
#' display=c(elecsoz,77:80)
#' heatmap_frag(frag=fragm3sp5s,elecsoz,time_window=c(-3,5),display=display)
heatmap_frag<-function(frag,elecsoz,time_window,option=NULL,title="PT01 seizure 1",display=display){
  
  titlepng=title
  fragdisplay=frag[display,]
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
#' @param scaling Scaling factor
#' @param displayChannels Channels to display
#'
#' @return plot raw signal
#' @export
#'
#' @examples
#' data("PT01Epochm30sp30s")
#' data("ElectrodesDataPT01")
#' displayChannels=which(ElectrodesDataPT01$insoz==TRUE)
#' visuiEEGdata(ieegts=PT01Epochm30sp30s,1000000, displayChannels = displayChannels)
visuiEEGdata<-function( ieegts, scaling, displayChannels){
  
  
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
#' @param quantilematrixsozsozc quantile matrix for the two groupd
#' @param time_window_ictal time window of processed ictal iEEG
#' @param subject_code subject code
#' @param j seizure number
#'
#' @return Quantile plot
#' @export
#'
#' @examples
#' data("fragstat")
#' time_window_ictal=c(-10,10)
#' plot_frag_quantile( quantilematrixsozsozc=fragstat[[1]], time_window_ictal=time_window_ictal,subject_code='pt01',j=1)
plot_frag_quantile<-function( quantilematrixsozsozc, time_window_ictal, subject_code,j){
 
  nw=ncol(quantilematrixsozsozc)
  stimes=c(1:nw)*(time_window_ictal[2]-time_window_ictal[1])/nw+time_window_ictal[1] 
  quantilesname<-c("SOZ(10th)","SOZ(20th)","SOZ(30th)","SOZ(40th)","SOZ(50th)",
                   "SOZ(60th)","SOZ(70th)","SOZ(80th)","SOZ(90th)","SOZ(100th)",
                   "SOZc(10th)","SOZc(20th)","SOZc(30th)","SOZc(40th)","SOZc(50th)",
                   "SOZc(60th)","SOZc(70th)","SOZc(80th)","SOZc(90th)","SOZc(100th)")
  quantileplot<- expand.grid(Time = stimes, Stats=quantilesname)
  quantileplot$Value <- c(t(quantilematrixsozsozc))
  
  titlepng=paste(subject_code,'Seizure',as.character(j),sep=" ")
  
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
#' @param cmeansoz time mean of soz group
#' @param cmeansozc time mean of non soz group
#' @param csdsoz time standard deviation of soz group
#' @param csdsozc time standard deviation of non soz group
#' @param time_window_ictal time window of processed ictal iEEG
#' @param subject_code subject code
#' @param j seizure number
#'
#' @return plot fragility distribution
#' @export
#'
#' @examples
#' data("fragstat")
#' time_window_ictal=c(-10,10)
#' plot_frag_distribution( cmeansoz=fragstat[[2]],cmeansozc=fragstat[[3]],csdsoz=fragstat[[4]],csdsozc=fragstat[[5]],time_window_ictal=time_window_ictal,subject_code='PT01',j=1)

plot_frag_distribution<-function( cmeansoz,cmeansozc,csdsoz,csdsozc,time_window_ictal,subject_code,j){
    
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
  
  titlepng=paste(subject_code,'Seizure',as.character(j),sep=" ")
  
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
