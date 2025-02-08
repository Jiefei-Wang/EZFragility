frag <- fg <- NULL
stat <- fstat <- NULL
nelec <- 10
ntime <- 100
elecsoz <- c(1, 2, 3)
ieegts <- pt01Fragility@ieegts[751:1250,]

test_that("calc_adj_frag", {
  set.seed(123)
  ieegts <- matrix(rnorm(ntime * nelec, -10, 10), ncol = nelec)
  t_window <- 20
  t_step <- 10
  frag <<- calc_adj_frag(
    ieegts = ieegts,
    t_window = t_window,
    t_step = t_step,
    nSearch = 2
  ) |> expect_no_error()
  expect_s4_class(frag, "Fragility")

  ## Test the show method
  print(frag) |> capture.output() |> expect_no_error()
})

test_that("GetAdjFrag", {
  capture.output({
    fg <<- GetAdjFrag(ieegts, 250, 125, 1e-4)
  }) |> expect_no_error()
  R2 <- fg$R2 |> expect_no_error()
  fg@R2 <- NULL
  print(fg) |> capture.output() |> expect_no_error()
  fg@R2 <- R2 |> expect_no_error()
})

test_that("frag_stat", {
  skip_if(!is(frag, "Fragility"))
  stat <<- frag_stat(frag = frag, elecsoz = elecsoz) |> expect_no_error()
  expect_s4_class(stat, "FragStat")
  ## Test the show method
  print(stat) |> capture.output() |> expect_no_error()
})

test_that("fragStat", {
  fstat <<- fragStat(fg, elecsoz) |> expect_no_error()
  qmat <- fstat$qmatrix |> expect_no_error()
  (fstat$qmatrix <- qmat) |> expect_no_error()
  fg@frag |> as.data.frame() |> fragStat(1:5) |> expect_error()
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
            if (length(nchar(n))) l[[n]] <- x else unNamed <- c(unNamed, x)
        }
        c(unNamed, l)
    }
}

i <- 77:84
iE <- 77:85
s <- colnames(ieegts)[i]
sE <- c(s, "Whatever")
soz <- 53:56

test_that("heatmap_frag", {
  dargs <- list(frag = fg, elecsoz = soz, time_window = c(-1, 2), title = "")
  vL <- def(dargs)
  do.call(heatmap_frag, dargs)  |> expect_no_error()
  do.call(heatmap_frag, vL(i))  |> expect_no_error()
  do.call(heatmap_frag, vL(iE)) |> expect_no_error()
  do.call(heatmap_frag, vL(s))  |> expect_no_error()
  do.call(heatmap_frag, vL(sE)) |> expect_no_error()
  do.call(heatmap_frag, vL(elecsoz = 85)) |> expect_no_error()
  do.call(heatmap_frag, vL(elecsoz = i)) |> expect_no_error()
  do.call(heatmap_frag, vL(elecsoz = iE)) |> expect_no_error()
  do.call(heatmap_frag, vL(elecsoz = s)) |> expect_no_error()
  do.call(heatmap_frag, vL(elecsoz = sE)) |> expect_no_error()
  do.call(heatmap_frag, vL(elecsoz = c("AD1", "AD2"))) |> expect_no_error()
})

test_that("visu_iEEG_data", {
  dargs <- list(
    ieegts = ieegts,
    time_window = c(-1, 1) / 2,
    title = ""
  )
  vL <- def(dargs)
  do.call(visu_iEEG_data, dargs) |> expect_no_error()
  do.call(visu_iEEG_data, vL(i)) |> expect_no_error()
  do.call(visu_iEEG_data, vL(iE)) |> expect_error()
  do.call(visu_iEEG_data, vL(s)) |> expect_no_error()
  do.call(visu_iEEG_data, vL(sE)) |> expect_no_error()
  do.call(visu_iEEG_data, vL(c("G1", "G2"))) |> expect_no_error()
})

test_that("plot_frag_distribution", {
  fstat |> plot_frag_distribution() |> expect_no_error()
})

test_that("plot_frag_quantile", {
  fstat |> plot_frag_quantile() |> expect_no_error()
})


## TODO:
# export(plot_frag_distribution)
# export(plot_frag_quantile)
# export(visuiEEGdata)
