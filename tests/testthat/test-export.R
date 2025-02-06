frag = fragObj <-  NULL
stat = fstat   <- NULL
nelec <- 10
ntime <- 100
elecsoz <- c(1, 2, 3)
sozindex <- attr(pt01Epoch, "sozindex")


test_that("calc_adj_frag", {
    set.seed(123)
    ieegts <- matrix(rnorm(ntime*nelec,-10, 10), ncol = nelec)
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
        fragObj <<- GetAdjFrag(matrix(rnorm(1e3), 100), 50, 25, 1e-4)
    }) |> expect_no_error()
    R2 <- fragObj$R2 |> expect_no_error()
    fragObj$R2 <- R2 |> expect_no_error()

})

test_that("frag_stat", {
    skip_if(!is(frag, "Fragility"))
    stat <<- frag_stat(frag = frag, elecsoz = elecsoz) |> expect_no_error()
    expect_s4_class(stat, "FragStat")
    ## Test the show method
    print(stat) |> capture.output() |> expect_no_error()

})

test_that("fragStat", {
    fstat <- fragStat(fragObj, elecsoz) |> expect_no_error()
    qmat <- fstat$qmatrix |> expect_no_error()
    (fstat$qmatrix <- qmat) |> expect_no_error()
    fragObj@frag |> as.data.frame() |> fragStat(1:5) |> expect_error()
})

test_that("heatmap_frag", {
    heatmap_frag(
        pt01Fragility,
        sozindex,
        c(-1,2),
        " ",
        c(sozindex, 77:80)
        ) |> expect_no_error()
})

test_that("heatmap_frag", {
    heatmap_frag(
        frag        = pt01Fragility,
        elecsoz     = sozindex,
        time_window = c(-1,2),
        title       = " ",
        display     = c(sozindex, 77:80)
    ) |> expect_no_error()
})

test_that("visu_iEEG_data", {
    visu_iEEG_data(
        ieegts      = pt01Epoch[9001:12000, ],
        time_window = c(-1, 2),
        title       = " ",
        display     = c(sozindex, 77:80)
    ) |> expect_no_error()
})


## TODO:
# export(plot_frag_distribution)
# export(plot_frag_quantile)
# export(visuiEEGdata)

