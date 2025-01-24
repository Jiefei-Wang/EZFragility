#' Visualization functions (raw signal, fragility matrix)
#' 
#' plot fragility heatmaps with electrodes marked as soz colored
#' 
#' @inheritParams frag_stat
#' @param elecsoz Integer or string. Vector soz electrodes (for good electrodes)
#' @param time_window Numeric Vector of length 2. The time window to display at the x-axis
#' @param title String. Figure title
#' @param display Integer or string. Vector electrodes to display
#'
#' @return Heatmap plot of the fragility matrix with soz electrodes in blue in the bottom
#'
#' @examples
#' # use integer index for display and soz electrodes
#'data("fragm3sp5s")
#'sozindex<-attr(fragm3sp5s,"sozindex")
#'time_window <- c(-3,5)
#'display <- c(sozindex,77:80)
#'fragplot<-heatmap_frag(frag=fragm3sp5s,elecsoz=sozindex,time_window <- time_window,title="PT01 seizure 1",display=display)
#'fragplot
#'
#' # use electrodes name for display and soz electrodes
#'data("fragm3sp5s")
#'soznames<-attr(fragm3sp5s,"soznames")
#'time_window <- c(-3,5)
#'display <- c(soznames,"MLT1","MLT2","MLT3","MLT4")
#'fragplot<-heatmap_frag(frag=fragm3sp5s,elecsoz=soznames,time_window <- time_window,title="PT01 seizure 1",display=display)
#'fragplot
#'
#' # save plot to file with ggplot2
#'data("fragm3sp5s")
#'sozindex<-attr(fragm3sp5s,"sozindex")
#'time_window <- c(-3,5)
#'display <- c(sozindex,77:80)
#'pathplot <- "~"
#'title <- "PT01sz1"
#'resfile <- paste(pathplot,'/FragilityHeatMap',title,'.png',sep="")
#'fragplot<-heatmap_frag(frag=fragm3sp5s,elecsoz=sozindex,time_window=time_window,title=title,display=display)
#'fragplot
#'ggplot2::ggsave(resfile)
#' 
#' @export
heatmap_frag<-function(frag,elecsoz,time_window = NULL,title="Patient name seizure number",display=NULL){
  titlepng<-title
  
  
  if(is.null(display)){
    display<-1:nrow(frag)
  }
  
  if (is(frag, "Fragility")) {
    frag <- frag$frag
  }

  
  elecname<-rownames(frag)
  elecind=c(1:nrow(frag))
  
  
  if(typeof(display)=="integer"){  
    
    displaytot<-1:nrow(frag)
    diffdisplaytot<-setdiff(display,displaytot)
    
    if(length(diffdisplaytot)!=0){
      listdisplaymissing<-paste(as.character(diffdisplaytot),collapse=" ")
      message<-paste("ERROR in display electrodes indices. Number(s) ",listdisplaymissing,"are out of electrode number limit")
      warning(message)
      display<-display[!display%in%diffdisplaytot]
      displaycor<-paste(as.character(display),collapse=" ")
      message<-paste("Keeping indices.",displaycor)
      warning(message)
      
    }  
    displayid<-display
    
  }else{
    
    diffdisplaytot<-setdiff(display,elecname)
    
    if(length(diffdisplaytot)!=0){
      listdisplaymissing<-paste(diffdisplaytot,collapse=" ")
      message<-paste("ERROR in display electrodes names. Name(s) ",listdisplaymissing,"are out of name list")
      warning(message)
      display<-display[!display%in%diffdisplaytot]
      displaycor<-paste(display,collapse=" ")
      message<-paste("Keeping names.",displaycor)
      warning(message)
      
    }  
    
    displayid<-which(elecname%in%display)
    
  }
  
  fragdisplay<-frag[displayid,]
  n_elec <- nrow(fragdisplay)
  electot<-c(1:n_elec)
  
  if(typeof(elecsoz)=="integer"){  
    
    diffelecind<-setdiff(elecsoz,elecind)
    
    if(length(diffelecind)!=0){
      listelecmissing<-paste(as.character(diffelecind),collapse=" ")
      message<-paste("ERROR in soz electrodes indices. Number(s) ",listelecmissing,"are out of electrode number limit")
      warning(message)
      elecsoz<-elecsoz[!elecsoz%in%diffelecind]
      listsozcor<-paste(as.character(elecsoz),collapse=" ")
      message<-paste("Keeping indices.",listsozcor)
      warning(message)
    }
    elecsozid<-elecsoz
    
  }else{

    diffsoztot<-setdiff(elecsoz,elecname)
    
    if(length(diffsoztot)!=0){
      listsozmissing<-paste(diffsoztot,collapse=" ")
      message<-paste("ERROR in soz electrodes names. Name(s) ",listsozmissing,"are out of name list")
      warning(message)
      elecsoz<-elecsoz[!elecsoz%in%diffsoztot]
      sozcor<-paste(elecsoz,collapse=" ")
      message<-paste("Keeping names.",sozcor)
      warning(message)
      
    }  
    
    elecsozid<-which(elecname%in%elecsoz)
  }

  elecsozd<-which(displayid%in%elecsozid)
  elecsozcd<-which(!displayid%in%elecsozid)
  elecsozsozc<-c(elecsozd,elecsozcd)

  elecnum <- rownames(fragdisplay)
  nw<- ncol(fragdisplay)
  colorelec<-elecnum
  nsoz<-length(elecsoz)
  colorelec[1:n_elec]<-"blue"
  nb<-n_elec-nsoz
  colorelec[1:nb]<-"black"

  elecsozsozc<-rev(elecsozsozc)
  fragord<-fragdisplay[elecsozsozc,]
  fragdf<-data.frame(fragord)
  
  if(is.null(time_window)){
    xlabel<-"Time Index"
    stimes<-seq_len(nw)
  }
  else{
    xlabel<-"Time (s)"
    stimes<-seq(time_window[1],time_window[2],length.out=nw)
  }

  colnames(fragdf)<-stimes
  rownames(fragdf)<-elecnum
  
  elecnum<-rev(elecnum)
  
  fragmap_data <- expand.grid(Time = stimes, Electrode = elecnum)
  fragmap_data$Value <- c(t(fragord))

  
  p<-ggplot2::ggplot(fragmap_data, ggplot2::aes(x = Time, y = Electrode, fill = Value)) +
    ggplot2::geom_tile() +
    ggplot2::ggtitle(as.character(titlepng)) +
    ggplot2::theme(plot.title=ggtext::element_markdown(hjust=0.5)) +
    ggplot2::labs(x = xlabel, y = "Electrode",size=2) +
    viridis::scale_fill_viridis(option = "turbo") +
    ggplot2::geom_vline(xintercept =0, 
                  color = "black", linetype = "dashed", size = 1) +
   ggplot2::theme_minimal() +
   ggplot2::theme(
      axis.text.y = ggplot2::element_text(size=6,colour=colorelec),     # Adjust depending on electrodes
    )
  
  return(p)

}

#' Visualization of ictal iEEG 
#'
#' @inheritParams heatmap_frag
#' @param ieegts Matrix or Fragility object. Either a matrix of iEEG time 
#' series x(t), with time points as rows and electrodes names as columns, 
#' or a Fragility object from \code{calc_adj_frag}
#' @param title String. Figure title
#' @param display Integer or string. Vector electrodes to display
#' @return plot raw signal
#'
#' @examples
#'data("pt01Epochm3sp5s")
#'sozindex<-attr(pt01Epochm3sp5s,"sozindex")
#'display <- c(sozindex,77:80)
#'time_window <- c(-3,5)
#'iEEGplot<-visu_iEEG_data(ieegts=pt01Epochm3sp5s,time_window=time_window,display=display)
#'iEEGplot
#' @export
visu_iEEG_data<-function(ieegts, time_window=NULL, title = "Patient name seizure number", display=NULL){
 
  titlepng<- title
  if(is.null(display)){
    display<-1:ncol(ieegts)
  }

  
  elecname<-colnames(ieegts)
  if(typeof(display)=="integer"){  

    displaytot<-1:nrow(ieegts)
    diffdisplaytot<-setdiff(display,displaytot)
    
    if(length(diffdisplaytot)!=0){
      listdisplaymissing<-paste(as.character(diffdisplaytot),collapse=" ")
      message<-paste("ERROR in display electrodes indices. Numbers ",listdisplaymissing,"are out of electrode number limit")
      warning(message)
      display<-display[!display%in%diffdisplaytot]
      displaycor<-paste(as.character(display),collapse=" ")
      message<-paste("Keeping indices.",displaycor)
      warning(message)
      
    }  
  
     displayid<-display
     
  }else{

    diffdisplaytot<-setdiff(display,elecname)
    
    if(length(diffdisplaytot)!=0){
      listdisplaymissing<-paste(diffdisplaytot,collapse=" ")
      message<-paste("ERROR in display electrodes names. Names ",listdisplaymissing,"are out of name list")
      warning(message)
      display<-display[!display%in%diffdisplaytot]
      displaycor<-paste(display,collapse=" ")
      message<-paste("Keeping names.",displaycor)
      warning(message)
      
    }  
        
    displayid<-which(elecname%in%display)
  }
  
  scaling <- 10^floor(log10(max(ieegts)))
  plotData<-ieegts[,displayid]/scaling
  gaps<-2
  displayNames<-colnames(ieegts)[displayid]
  n_elec<-length(displayid)
  nt<-nrow(plotData)
  if(is.null(time_window)){
    xlabel<-"Time Index"
    stimes<-seq_len(nt)
  }
  else{
    xlabel<-"Time (s)"
    stimes<-seq(time_window[1],time_window[2],length.out=nt)
  }


  for(i in 1:ncol(plotData)){
     plotData[, i] <- (plotData[, i]- mean(plotData[, i]))+
       (ncol(plotData)-i)*gaps
  }
  plotData<-data.frame(plotData)
  breakplot<-(c(1:n_elec)-1)*gaps
 
  p<-ggplot2::ggplot(data=plotData,ggplot2::aes(x=stimes,y=plotData))+
  ggplot2::ggtitle(titlepng)+
  ggplot2::labs(x = xlabel, y = "Electrode",size=2)+ 
    ggplot2::geom_vline(xintercept =0, 
                  color = "black", linetype = "dashed", size = 1)
  
  for(i in 1:n_elec){
      p<-p+ggplot2::geom_line(ggplot2::aes_string(y=names(plotData)[i]))
  }
  displayNames<-rev(displayNames)
  p<-p+ggplot2::scale_y_continuous(labels=displayNames,breaks=breakplot)
  
  return(p)

}

#' Plot Fragility time quantiles for two electrodes group marked as soz non marked as soz
#'
#' @inheritParams heatmap_frag
#' @param qmatrix Matrix or FragStat object, either a quantile matrix 
#' for the two groups or a FragStat object from \code{frag_stat}
#' @param title String. Figure title
#'
#' @return Quantile plot
#' @export
#'
#' @examples
#' time_window <- c(-3,5)
#'data("fragm3sp5s")
#'sozindex<-attr(fragm3sp5s,"sozindex")
#'# compute fragility statistics evolution with time (mean and standard deviation) for soz and
#'# non soz groups
#'fragstat <- frag_stat(frag=fragm3sp5s, elecsoz=sozindex)
#'qfragplot<-plot_frag_quantile(stat=fragstat, time_window=time_window)
#'qfragplot
plot_frag_quantile<-function(stat, time_window = NULL,title="Fragility Quantiles over time"){
  
  if(!is.null(stat)){
    qmatrix<-stat$qmatrix
  }
  nw <- ncol(qmatrix)
  if(is.null(time_window)){
    xlabel<-"Time Index"
    stimes<-seq_len(nw)
  }
  else{
    xlabel<-"Time (s)"
    stimes<-seq(time_window[1],time_window[2],length.out=nw)
  }


  quantilesname<-c(paste0("SOZ(", seq(10,100,by=10),")"),paste0("SOZc(", seq(10,100,by=10),")"))
  quantileplot<- expand.grid(Time = stimes, Stats=quantilesname)
  quantileplot$Value <- c(t(qmatrix))
  
  titlepng <- title
  
  p<-ggplot2::ggplot(quantileplot, ggplot2::aes(x = Time, y = Stats, fill = Value)) +
    ggplot2::geom_tile() +
    ggplot2::ggtitle(titlepng)+
    ggplot2::labs(x = xlabel, y = "Quantiles",size=2) +
    viridis::scale_fill_viridis(option = "turbo") +  #
    
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size=4),     # Adjust depending on electrodes
    )

  if(!is.null(time_window)){
    p <- p + ggplot2::geom_vline(xintercept =0, 
                                 color = "black", linetype = "dashed", size = 1)}

  return(p)
}


#'  Plot Fragility time distribution for two electrodes group marked and non-marked as soz
#'
#' @inheritParams heatmap_frag
#' @param stat FragStat object, a FragStat object from \code{frag_stat}, if specified, the arguments cmeansoz, cmeansozc, csdsoz, csdsozc will be ignored
#' @param cmeansoz Numeric Vector. mean soz group in function of time window
#' @param cmeansozc Numeric Vector. mean non soz group in function of time window
#' @param csdsoz Numeric Vector. standard deviation soz group in function of time window
#' @param csdsozc Numeric Vector. standard deviation non soz group in function of time window
#' @param title String. Figure title
#'
#' @return plot fragility distribution
#' @export
#'
#' @examples
#'time_window <- c(-3,5)
#'data("fragm3sp5s")
#'sozindex<-attr(fragm3sp5s,"sozindex")
#'# compute fragility statistics evolution with time (mean and standard deviation) for soz and
#'# non soz groups
#'fragstat <- frag_stat(frag=fragm3sp5s, elecsoz=sozindex)
#'# plot the statistical results
#'pfragstat<-plot_frag_distribution(stat=fragstat,time_window=time_window)
#'pfragstat
plot_frag_distribution<-function(
  stat = NULL, 
  time_window = NULL,
  cmeansoz = NULL, cmeansozc = NULL, csdsoz = NULL, csdsozc = NULL,
  title='Average Fragility over time'){
  if(!is.null(stat)){
    cmeansoz <- stat$cmeansoz
    cmeansozc <- stat$cmeansozc
    csdsoz <- stat$csdsoz
    csdsozc <- stat$csdsozc
  }
  
  nw <- length(cmeansoz)
  if(is.null(time_window)){
    xlabel<-"Time Index"
    stimes<-seq_len(nw)
  }
  else{
    xlabel<-"Time (s)"
    stimes<-seq(time_window[1],time_window[2],length.out=nw)
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
  p<-ggplot2::ggplot(plotmeanstd, ggplot2::aes(x=times, y=cmeansoz))+ 
   ggplot2::xlab(xlabel)+
   ggplot2::ylab('Fragility')+
   ggplot2::ggtitle(titlepng)+
   ggplot2::geom_line(ggplot2::aes(y = meansoz,color="SOZ +/- sem"))+  
   ggplot2::geom_line(ggplot2::aes(y = sozsdp),color='red',linetype="dotted")+  
   ggplot2::geom_line(ggplot2::aes(y = sozsdm),color='red',linetype="dotted")+ 
   ggplot2::geom_line(ggplot2::aes(y = meansozc,color="SOZc +/- sem"))+
   ggplot2::geom_line(ggplot2::aes(y = sozcsdp),color='black',linetype="dotted")+  
   ggplot2::geom_line(ggplot2::aes(y = sozcsdm),color='black',linetype="dotted")+
   ggplot2::geom_ribbon(ggplot2::aes(ymin=sozsdm,ymax=sozsdp), fill="red",alpha=0.5)+
   ggplot2::geom_ribbon(ggplot2::aes(ymin=sozcsdm,ymax=sozcsdp), fill="black",alpha=0.5)+  
   ggplot2::scale_color_manual(name="Electrode groups",values = c(colors)) 
    
  ## add vertical line at time 0 if time_window is specified
  if(!is.null(time_window)){
    p <- p + ggplot2::geom_vline(xintercept =0, 
                  color = "black", linetype = "dashed", size = 1)}
  return(p)
}
