context("mcolnames")

test_that("DataFrame", {
    expect_identical(
        object = mcols(df),
        expected = NULL
    )
})

test_that("GRanges", {
    expect_identical(
        object = mcolnames(gr),
        expected = c("geneID", "geneName")
    )
})

test_that("SummarizedExperiment", {
    expect_true("geneName" %in% mcolnames(rse))
})



context("mcolnames<-")

test_that("Don't allow length mismatch", {
    expect_error(
        object = mcolnames(df) <- "xxx",
        regexp = "areSameLength"
    )
})
