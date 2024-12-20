#' Visualization
#' 
#' plot fragility heatmaps with electrodes marked as soz colored
#'
#' @param frag fragility matrix results
#' @param goodelec 
#' @param sozelec 
#' @param ieegts Numeric. A matrix of iEEG time series x(t), 
#' with time points as rows and electrodes names as columns
#' @param time_window Fragility heatmap time window around seizure onset
#' @param subject_code patient name
#' @param j seizure number
#'
#' @return a heatmap plot of the fragility matrix with soz electrodes in blue in the bottom
#' @export
#'
#' @examples
#' data("fragm3sp5s")
#' data("fragrankm3sp5s")
#' data("goodelec")
#' data("pt01Epochm3sp5s")
#' data("sozelec")
#' time_window=c(-3:5)
#' heatmap_frag(frag=fragm3sp5s,goodelec=goodelec,sozelec=sozelec,ieegts=pt01Epochm3sp5s,time_window=c(-3,5),subject_code='pt01',j=1)
heatmap_frag<-function(frag,goodelec,sozelec,ieegts,time_window,option=NULL,subject_code,j){
  
  titlepng=paste(subject_code,'Seizure',as.character(j),sep=" ")
  insoz=goodelec%in%sozelec
  elecsoz=which(insoz==TRUE)
  elecsozc=which(insoz==FALSE)
  elecsozsozc=c(elecsoz,elecsozc)
  
  elecnum <- colnames(ieegts)[elecsozsozc]
  n_elec <- ncol(ieegts)
  nw<- ncol(frag)
  colorelec<-elecnum
  nsoz=length(elecsoz)
  colorelec[1:n_elec]="black"
  colorelec[1:nsoz]="blue"

  fragord<-frag[elecsozsozc,]
  fragdf<-data.frame(fragord)
  stimes=c(1:nw)*(time_window[2]-time_window[1])/nw+time_window[1]
  colnames(fragdf)<-stimes
  rownames(fragdf)<-elecnum
  
  fragmap_data <- expand.grid(Time = stimes, Electrode = elecnum)
  fragmap_data$Value <- c(t(fragord))
  
  ggplot2::ggplot(fragmap_data, ggplot2::aes(x = Time, y = Electrode, fill = Value)) +
    ggplot2::geom_tile() +
    ggplot2::ggtitle(titlepng)+
    ggplot2::labs(x = "Time (s)", y = "Electrode") +
    viridis::scale_fill_viridis(option = "turbo") +  #
    
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size=5,colour=colorelec),     # Adjust depending on electrodes
    )
  
}