#' Visualization
#' 
#' plot fragility heatmaps with electrodes marked as soz colored
#'
#' @param frag fragility matrix results
#' @param goodelec 
#' @param sozelec 
#' @param ieegts Numeric. A matrix of iEEG time series x(t), 
#' with time points as rows and electrodes names as columns
#'
#' @return
#' @export
#'
#' @examples
#' heatmap_frag(frag=frag,goodelec=goodelec,sozelec=sozelec,ieegts=ptEpochm1sp2s,subject_code='pt01',j=1)
heatmap_frag<-function(frag,goodelec,sozelec,ieegts,option=NULL,plotelec,subject_code,j){
  
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
  stimes=c(1:nw)*3/nw-1
  colnames(fragdf)<-stimes
  rownames(fragdf)<-elecnum
  
  fragmap_data <- expand.grid(Time = stimes, Electrode = elecnum)
  fragmap_data$Value <- c(t(fragord))
  
 # library(ggplot2) 
 # library(viridis)
  ggplot(fragmap_data, aes(x = Time, y = Electrode, fill = Value)) +
    geom_tile() +
    ggtitle(titlepng)+
    labs(x = "Time (s)", y = "Electrode") +
    scale_fill_viridis(option = "turbo") +  #
    
    theme_minimal() +
    theme(
      axis.text.y = element_text(size=5,colour=colorelec),     # Adjust depending on electrodes
    )
  
}