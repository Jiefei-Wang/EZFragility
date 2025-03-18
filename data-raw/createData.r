
## load pt01epochdata.mat
## Patient PT01 from the Fragility data set

library(R.matlab)
library(readxl)
data <- readMat('data-raw/pt01epochdata.mat')
pt01EpochRaw <- data$a

## add channel names to the rows
goodChannels <- c(1:4,7:36,42:43,46:69,72:95)
sozChannels<-c(33:34,62:69)
channelNames <- read_excel('data-raw/Pt01ictalRun01EcoGChannels.xls')
rownames(pt01EpochRaw) <- channelNames$name[goodChannels]
sozIndex<-which(goodChannels%in%sozChannels==TRUE)
sozNames<-channelNames$name[sozChannels]

## Add time stamps to the columns
times <- seq(-10, 10, length.out=ncol(pt01EpochRaw))
times_with_sign <- ifelse(times >= 0, paste0("+", times), as.character(times))
colnames(pt01EpochRaw)<-times_with_sign


pt01EcoG<-pt01EpochRaw[,9001:12000]
attr(pt01EcoG, "sozIndex") <- sozIndex
attr(pt01EcoG, "sozNames") <- sozNames
usethis::use_data(pt01EcoG, overwrite = TRUE)


## load fragility matrix
# library(doSNOW)
# library(EZFragility)
# cl <- makeCluster(8, type = "SOCK")
# registerDoSNOW(cl)
t_window <- 250
t_step <- 125
pt01Frag <- calcAdjFrag(epoch = pt01EcoG, window = t_window, step = t_step, parallel = FALSE, progress = TRUE)
# stopCluster(cl)
usethis::use_data(pt01Frag, overwrite = TRUE)


