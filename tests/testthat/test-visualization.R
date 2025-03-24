data(pt01EcoG)
data(pt01Frag)
displayIndex <- 77:84
displayIndexOutOfRange <- 77:85
str <- rownames(pt01EcoG)[displayIndex]
strError <- c(str, "Whatever")

test_that("visuIEEGData", { 
  visuIEEGData(pt01EcoG) |> expect_no_error()
})

test_that("plotFragHeatmap", {
  plotFragHeatmap(pt01Frag, displayIndex)           |> expect_no_error()
  plotFragHeatmap(pt01Frag, displayIndexOutOfRange) |> expect_warning()
  plotFragHeatmap(pt01Frag, str)                    |> expect_no_error()
  plotFragHeatmap(pt01Frag, strError)               |> expect_warning()
})

test_that("plotFragDistribution", {
  plotFragDistribution(pt01Frag, displayIndex)           |> expect_no_error()
  plotFragDistribution(pt01Frag, displayIndexOutOfRange) |> expect_warning()
  plotFragDistribution(pt01Frag, str)                    |> expect_no_error()
  plotFragDistribution(pt01Frag, strError)               |> expect_warning()
})

test_that("plotFragQuantile", {
  plotFragQuantile(pt01Frag, displayIndex)           |> expect_no_error()
  plotFragQuantile(pt01Frag, displayIndexOutOfRange) |> expect_warning()
  plotFragQuantile(pt01Frag, str)                    |> expect_no_error()
  plotFragQuantile(pt01Frag, strError)               |> expect_warning()
})
