set.seed(1)
data <- matrix(rnorm(800), ncol = 40)
window <- 10
step <- 5
lambda <- 0.1
testFrag <- calcAdjFrag(epoch = data, window = window, step = step, lambda = lambda)
usethis::use_data(testFrag, overwrite = TRUE, internal = TRUE)
