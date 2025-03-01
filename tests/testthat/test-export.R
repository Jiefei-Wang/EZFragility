frag <- NULL
stat <-  NULL
nelec <- 10
ntime <- 100
elecSoz <- c(1, 2, 3)

test_that("calcAdjFrag", {
  set.seed(123)
  ieegts <- matrix(rnorm(ntime * nelec, -10, 10), ncol = nelec)
  window <- 20
  step <- 10
  frag <<- calcAdjFrag(
    ieegts = ieegts,
    window = window,
    step = step,
    nSearch = 2
  ) |> expect_no_error()
  expect_s4_class(frag, "Fragility")

  ## Test the show method
  print(frag) |> capture.output() |> expect_no_error()
})

test_that("fragStat", {
  skip_if(!is(frag, "Fragility"))
  stat <<- fragStat(frag = frag, sozID = elecSoz) |> expect_no_error()
  expect_s4_class(stat, "FragStat")
  ## Test the show method
  print(stat) |> capture.output() |> expect_no_error()
})


test_that("S4 methods", {
  frag@lambdas <- NULL
  print(frag) |> capture.output() |> expect_no_error()
  
  ## Test assignment
  frag |> is("Fragility") |> expect_equal(TRUE)
  frag$lambdas <- NULL
  frag |> is("Fragility") |> expect_equal(TRUE)
})

# Visualization ----------------------------------------------------------------
def <- \(x) {
    l <- x
    \(...) {
        inp <- list(...)
        unNamed <- list()
        for (i in seq_along(inp)) {
            n <- names(inp[i])
            x <- if (is.symbol(inp[[i]])) eval(inp[[i]]) else inp[[i]]
            if (is.atomic(x)) x <- list(x)
            if (length(nchar(n))) l[n] <- x
            else unNamed <- c(unNamed, if (is.atomic(x)) list(x) else x)
        }
        c(unNamed, l)
    }
}

data(pt01Epochm1sp2s)
data(pt01Fragm1sp2s)
fg <- pt01Fragm1sp2s
int <- 77:84
intError <- 77:85
str <-  colnames(pt01Epochm1sp2s)[int]
strError <- c(str, "Whatever")
soz <- 53:56

test_that("valid_soz", {
  mat <- pt01Epochm1sp2s
  valid_soz(mat, int, str) |> expect_no_error()
  valid_soz(mat, intError, strError) |> expect_warning() |> expect_warning()
})

test_that("heatmapFrag", {
  dargs <- list(frag = fg, sozID = soz, timeRange = c(-1, 2), title = "")
  vL <- def(dargs)
  do.call(heatmapFrag, dargs) |> expect_no_error()
  do.call(heatmapFrag, vL(timeRange = NULL)) |> expect_no_error()
  
  do.call(heatmapFrag, vL(int)) |> expect_no_error()
  do.call(heatmapFrag, vL(intError)) |> expect_warning() |> expect_warning()
  
  do.call(heatmapFrag, vL(str))      |> expect_no_error()
  do.call(heatmapFrag, vL(strError)) |> expect_warning() |> expect_warning()
  
  do.call(heatmapFrag, vL(sozID = int))      |> expect_no_error()
  do.call(heatmapFrag, vL(sozID = intError)) |> expect_warning()  |> expect_warning()
  
  do.call(heatmapFrag, vL(sozID = str))      |> expect_no_error()
  do.call(heatmapFrag, vL(sozID = strError)) |> expect_warning()  |> expect_warning()
})

test_that("visuIEEGData", {
  dargs <- list( ieegts = pt01Epochm1sp2s, timeRange = c(-1, 2), title = "" )
  vL <- def(dargs)
  do.call(visuIEEGData, dargs)        |> expect_no_error()
  do.call(visuIEEGData, vL(int))      |> expect_no_error()
  do.call(visuIEEGData, vL(intError)) |> expect_warning() |> expect_warning()
  do.call(visuIEEGData, vL(str))      |> expect_no_error()
  do.call(visuIEEGData, vL(strError)) |> expect_warning() |> expect_warning()
  do.call(visuIEEGData, vL(timeRange = NULL)) |> expect_no_error()
})

test_that("plotFragDistribution", {
  stat |> plotFragDistribution()         |> expect_no_error()
  stat |> plotFragDistribution(c(-1, 2)) |> expect_no_error()
})

test_that("plotFragQuantile", {
  stat |> plotFragQuantile()             |> expect_no_error()
  stat |> plotFragDistribution(c(-1, 2)) |> expect_no_error()
})
