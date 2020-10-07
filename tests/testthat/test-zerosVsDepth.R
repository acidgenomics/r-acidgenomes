context("zerosVsDepth")

test_that("SummarizedExperiment", {
    x <- zerosVsDepth(rse)
    expect_s4_class(x, "DataFrame")
    expect_identical(
        object = round(mean(x[["dropout"]]), digits = 3L),
        expected = 0.106
    )
    expect_identical(
        object = round(mean(x[["depth"]]), digits = 2L),
        expected = 25305.67
    )
})

test_that("SingleCellExperiment", {
    x <- zerosVsDepth(sce)
    expect_s4_class(x, "DataFrame")
    expect_identical(
        object = round(mean(x[["dropout"]]), digits = 2L),
        expected = 0.48
    )
    expect_identical(
        object = round(mean(x[["depth"]]), digits = 2L),
        expected = 57634.36
    )
})
