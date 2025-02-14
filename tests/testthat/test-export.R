frag <- NULL
stat <-  NULL
nelec <- 10
ntime <- 100
elecsoz <- c(1, 2, 3)

test_that("calc_adj_frag", {
  set.seed(123)
  ieegts <- matrix(rnorm(ntime * nelec, -10, 10), ncol = nelec)
  t_window <- 20
  t_step <- 10
  frag <<- calc_adj_frag(
    ieegts = ieegts,
    t_window = t_window,
    t_step = t_step,
    n_search = 2
  ) |> expect_no_error()
  expect_s4_class(frag, "Fragility")

  ## Test the show method
  print(frag) |> capture.output() |> expect_no_error()
})

test_that("frag_stat", {
  skip_if(!is(frag, "Fragility"))
  stat <<- frag_stat(frag = frag, elecsoz = elecsoz) |> expect_no_error()
  expect_s4_class(stat, "FragStat")
  ## Test the show method
  print(stat) |> capture.output() |> expect_no_error()
})


test_that("S4 methods", {
  frag@lambdas <- NULL
  print(frag) |> capture.output() |> expect_no_error()
  
  # Should this be allowed???
  frag@R2 <- matrix(LETTERS[1:24], 6)
  print(frag) |> capture.output() |> expect_no_error()
  stat1 <- stat
  stat1@csdsozc <- stat1$csdsozc[1:2]
  print(stat1) |> capture.output() |> expect_no_error()

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

data(pt01Frag)
fg <- pt01Frag
int <- 77:84
intError <- 77:85
str <-  colnames(fg@ieegts)[int]
strError <- c(str, "Whatever")
soz <- 53:56

test_that("valid_soz", {
  mat <- fg@ieegts
  valid_soz(mat, int, str) |> expect_no_error()
  valid_soz(mat, intError, strError) |> expect_warning() |> expect_warning()
})

test_that("heatmap_frag", {
  dargs <- list(frag = fg, elecsoz = soz, time_window = c(-1, 2), title = "")
  vL <- def(dargs)
  do.call(heatmap_frag, dargs) |> expect_warning()
  do.call(heatmap_frag, vL(time_window = NULL)) |> expect_warning()
  
  do.call(heatmap_frag, vL(int)) |> expect_warning()
    do.call(heatmap_frag, vL(intError)) |> expect_warning() |> 
    expect_warning() |>
    expect_warning()
  
  do.call(heatmap_frag, vL(str))      |> expect_warning()
  do.call(heatmap_frag, vL(strError)) |> expect_warning() |> 
    expect_warning() |> 
    expect_warning()
  
  do.call(heatmap_frag, vL(elecsoz = int))      |> expect_warning()
  do.call(heatmap_frag, vL(elecsoz = intError)) |> expect_warning() |> 
    expect_warning() |> 
    expect_warning()
  
  do.call(heatmap_frag, vL(elecsoz = str))      |> expect_warning()
    do.call(heatmap_frag, vL(elecsoz = strError)) |> expect_warning() |> 
    expect_warning() |> 
    expect_warning()
})

test_that("visu_iEEG_data", {
  dargs <- list( ieegts = fg$ieegts, time_window = c(-1, 2), title = "" )
  vL <- def(dargs)
  do.call(visu_iEEG_data, dargs)        |> expect_no_error()
  do.call(visu_iEEG_data, vL(int))      |> expect_no_error()
  do.call(visu_iEEG_data, vL(intError)) |> expect_error()
  do.call(visu_iEEG_data, vL(str))      |> expect_no_error()
  do.call(visu_iEEG_data, vL(strError)) |> expect_warning() |> expect_warning()
  do.call(visu_iEEG_data, vL(time_window = NULL)) |> expect_no_error()
})

test_that("plot_frag_distribution", {
  stat |> plot_frag_distribution()         |> expect_no_error()
  stat |> plot_frag_distribution(c(-1, 2)) |> expect_no_error()
})

test_that("plot_frag_quantile", {
  stat |> plot_frag_quantile()             |> expect_no_error()
  stat |> plot_frag_distribution(c(-1, 2)) |> expect_no_error()
})
