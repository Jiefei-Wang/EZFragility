frag <- NULL
stat <- NULL
nelec <- 10
ntime <- 100
elecsoz <- c(1, 2, 3)


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


test_that("frag_stat", {
    skip_if(!is(frag, "Fragility"))
    stat <<- frag_stat(frag = frag, elecsoz = elecsoz) |> expect_no_error()
    expect_s4_class(stat, "FragStat")

    ## Test the show method
    print(stat) |> capture.output() |> expect_no_error()
})

test_that("heatmap_frag", {
    skip_if(!is(frag, "Fragility"))
    ## Why this does not work?
    heatmap_frag(frag = frag, elecsoz = elecsoz) |> expect_no_error()
})


## TODO:
# export(plot_frag_distribution)
# export(plot_frag_quantile)
# export(visuiEEGdata)

