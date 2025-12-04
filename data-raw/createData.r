
## Patient PT01 from the Fragility data set

library(Epoch)
library(gsignal)

# download data
dl<-EpochDownloader(progress = FALSE)
names(dl)

#The preprocessed voltage data from patient pt01 seizure 1 can be loaded by
pt01sz1 <- dl$FragilityData_subpt01_1

pt01EcoG <- pt01EcoG |>
  crop(start=-1, end=2)
pt01EcoG

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




