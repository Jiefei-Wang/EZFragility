test_that("calc_adj_frag works", {
  set.seed(123)
  mat <- matrix(rnorm(1000,-10, 10), ncol = 10)
  calc_adj_frag(mat, 20, 10, 0.1) |> expect_no_error()
  set.seed(123)
  calc_adj_frag(mat, 20, 10) |> expect_no_error()
  ridge(mat[1:4, ], mat[1:4 + 1, ], 0.1, FALSE, NULL) |>
    fragilityRow() |>
    expect_no_error()
  ridge(mat[1:4, ], mat[1:4 + 1, -1], 0.1, FALSE, NULL) |> expect_error()
  ridgesearchlambdadichomotomy(mat[1:4, ], mat[1:4 + 1, -1], FALSE, NULL) |>
    expect_error()
  ridgesearchlambdadichomotomy(mat[1:4, ], mat[1:4 + 1, ], TRUE, NULL) |>
    expect_no_error()
})

test_that("frag_quantile works", {
  # frag_quantile(repository, f, t_window, t_step, soz, sozc) |> expect_no_error()
  # mean_f_calc(repository, f, soz, sozc) |> expect_no_error()
  # threshold_buckets(f, c(0.1, 0.5, 1)) |> expect_no_error()
  # threshold_fragility(repository, adj_frag_info, t_step, 0.1, 1) |>
  #   expect_no_error()
})
