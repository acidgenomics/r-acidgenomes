context("geneNames")

genes <- c("TSPAN6", "TNMD", "DPM1", "SCYL3", "C1orf112")

test_that("GRanges", {
    expect_identical(
        object = geneNames(gr),
        expected = genes
    )
})

test_that("SummarizedExperiment", {
    expect_identical(
        object = head(geneNames(rse), n = 5L),
        expected = genes
    )
})
