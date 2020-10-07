context("autopadZeros")

test_that("SummarizedExperiment", {
    x <- autopadZeros(
        object = rse,
        rownames = TRUE,
        colnames = TRUE,
        sort = TRUE
    )
    expect_s4_class(x, "RangedSummarizedExperiment")
})

test_that("SingleCellExperiment", {
    x <- autopadZeros(
        object = sce,
        rownames = TRUE,
        colnames = TRUE,
        sort = TRUE
    )
    expect_s4_class(x, "SingleCellExperiment")
})
