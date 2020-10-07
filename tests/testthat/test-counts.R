context("counts")

test_that("SummarizedExperiment", {
    object <- counts(rse)
    expect_is(object, "matrix")
})

test_that("SE assignment", {
    counts <- counts(rse)
    counts <- log2(counts + 1L)
    counts(rse) <- counts
    expect_s4_class(rse, "RangedSummarizedExperiment")
})

test_that("SE assignment failure", {
    expect_error(counts(rse) <- matrix())
})
