context("makeSingleCellExperiment")

assays <- assays(sce)

test_that("Minimal mode", {
    x <- makeSingleCellExperiment(assays = assays)
    expect_s4_class(x, "SingleCellExperiment")
})
