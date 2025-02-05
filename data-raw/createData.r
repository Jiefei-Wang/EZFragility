## load pt01epochdata.mat

library(R.matlab)
library(readxl)
data <- readMat('pt01epochdata.mat')
ptEpoch <- data$a

## add channel names to the rows
goodChannels <- c(1:4,7:36,42:43,46:69,72:95)
channelNames <- read_excel('data-raw/Pt01ictalRun01EcoGChannels.xls')
rownames(ptEpoch) <- channelNames$name[goodChannels]

## Add time stamps to the columns
times <- seq(-10, 10, length.out=ncol(ptEpoch))
times_with_sign <- ifelse(times >= 0, paste0("+", times), as.character(times))
colnames(ptEpoch) <- paste0('t', times_with_sign)

ptEpoch <- t(ptEpoch)

usethis::use_data(ptEpoch, overwrite = TRUE)

