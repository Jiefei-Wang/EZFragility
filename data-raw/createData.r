
## load pt01epochdata.mat
## Patient PT01 from the Fragility data set

library(R.matlab)
library(readxl)
data <- readMat('data-raw/pt01epochdata.mat')
pt01Epochm1sp2s <- data$a

## add channel names to the rows
goodChannels <- c(1:4,7:36,42:43,46:69,72:95)
sozChannels<-c(33:34,62:69)
channelNames <- read_excel('data-raw/Pt01ictalRun01EcoGChannels.xls')
rownames(pt01Epochm1sp2s) <- channelNames$name[goodChannels]
sozIndex<-which(goodChannels%in%sozChannels==TRUE)
sozNames<-channelNames$name[sozChannels]

## Add time stamps to the columns
times <- seq(-10, 10, length.out=ncol(pt01Epochm1sp2s))
times_with_sign <- ifelse(times >= 0, paste0("+", times), as.character(times))
colnames(pt01Epochm1sp2s)<-times_with_sign

pt01Epochm1sp2s <- t(pt01Epochm1sp2s)
pt01Epochm1sp2s<-pt01Epochm1sp2s[9001:12000,]
attr(pt01Epochm1sp2s, "sozIndex") <- sozIndex
attr(pt01Epochm1sp2s, "sozNames") <- sozNames
usethis::use_data(pt01Epochm1sp2s, overwrite = TRUE)


## load fragility matrix
t_window <- 250
t_step <- 125
lambda <- NULL
pt01Fragm1sp2s <- calcAdjFrag(ieegts = pt01Epochm1sp2s, window = t_window, step = t_step, lambda = lambda)
usethis::use_data(pt01Fragm1sp2s, overwrite = TRUE)


