data("pt01EcoG")
epoch <- Epoch(pt01EcoG)

test_that("truncateTime", {
  truncateTime(epoch, 0, 1) |> expect_no_error()
})
