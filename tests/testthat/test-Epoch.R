data("pt01EcoG")

test_that("Epoch", {
  epoch <- Epoch(pt01EcoG) |> expect_no_error()
  epoch |> expect_s4_class("Epoch")
})
