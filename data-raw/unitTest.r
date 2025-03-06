set.seed(1)
data <- matrix(rnorm(800), nrow = 40)
t_window <- 10
t_step <- 5
lambda <- 0.1
testFrag <- calcAdjFrag(ieegts = data, window = t_window, step = t_step, lambda = lambda)
usethis::use_data(testFrag, overwrite = TRUE, internal = TRUE)
