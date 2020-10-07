context("topCellsPerSample")

test_that("SingleCellExperiment", {
    x <- topCellsPerSample(sce, n = 5L)
    expect_is(x, "list")
    expect_identical(
        object = names(x),
        expected = paste0("sample", seq(2L))
    )
})
