test_that("calc_adj_frag works", {
  set.seed(123)
  mat <- matrix(rnorm(1000,-10, 10), ncol = 10)
  calc_adj_frag(mat, 20, 10, 0.1) |> expect_no_error()
  set.seed(123)
  calc_adj_frag(mat, 20, 10) |> expect_no_error()
  ridge(mat[1:4, ], mat[1:4 + 1, ], 0.1, FALSE) |> expect_no_error()
  A= matrix(
    c(-1,3,-3,-1), 
    nrow = 2,             
    ncol = 2,             
    byrow = TRUE          
  )
  fragilityRow(A) |> expect_no_error()
  ridgesearchlambdadichomotomy(mat[1:4, ], mat[1:4 + 1, -1], FALSE, NULL) |>
   expect_error()
  ridgesearchlambdadichomotomy(mat[1:4, ], mat[1:4 + 1, ], FALSE) |>
   expect_no_error()
})

test_that("frag_stat works", {
  set.seed(123)
  mat <- matrix(rnorm(1000,-10, 10), ncol = 100)
  frag_stat(mat,c(1:5))|>
    expect_no_error()
})

