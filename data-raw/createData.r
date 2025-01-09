## load pt01epochdata.mat
## Patient PT01 from the Fragility data set

library(R.matlab)
library(readxl)
data <- readMat('data-raw/pt01epochdata.mat')
pt01Epoch <- data$a

## add channel names to the rows
goodChannels <- c(1:4,7:36,42:43,46:69,72:95)
sozChannels<-c(33:34,62:69)
channelNames <- read_excel('data-raw/Pt01ictalRun01EcoGChannels.xls')
rownames(pt01Epoch) <- channelNames$name[goodChannels]
sozindex<-which(goodChannels%in%sozChannels==TRUE)
soznames<-channelNames$name[sozChannels]

## Add time stamps to the columns
times <- seq(-10, 10, length.out=ncol(pt01Epoch))
times_with_sign <- ifelse(times >= 0, paste0("+", times), as.character(times))
# why this?
#colnames(ptEpoch) <- paste0('t', times_with_sign)
colnames(pt01Epoch)<-times_with_sign

pt01Epoch <- t(pt01Epoch)
attr(pt01Epoch, "sozindex") <- sozindex
attr(pt01Epoch, "soznames") <- soznames
usethis::use_data(pt01Epoch, overwrite = TRUE)

pt01Epochm1sp2s<-pt01Epoch[9001:12000,]
attr(pt01Epochm1sp2s, "sozindex") <- sozindex
attr(pt01Epochm1sp2s, "soznames") <- soznames
usethis::use_data(pt01Epochm1sp2s, overwrite = TRUE)

pt01Epochm3sp5s<-pt01Epoch[7001:15000,]
attr(pt01Epochm3sp5s, "sozindex") <- sozindex
attr(pt01Epochm3sp5s, "soznames") <- soznames
usethis::use_data(pt01Epochm3sp5s, overwrite = TRUE)

