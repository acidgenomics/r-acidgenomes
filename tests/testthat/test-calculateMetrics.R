context("calculateMetrics")

test_that("matrix-like", {
    for (object in list(
        Matrix = counts(sce)
        ## nolint start
        ## > DelayedArray = DelayedArray(counts(sce))
        ## nolint end
    )) {
        object <- calculateMetrics(object, prefilter = TRUE)
        expect_s4_class(object, "DataFrame")
        ## Check for run-length encoding.
        expect_true(all(bapply(object, is, class2 = "Rle")))
        ## Check that expected values match.
        ## Parameterized unit test is changing `nCount` to NA.
        x <- object[1L, c("nCount", "nFeature", "nMito"), drop = TRUE]
        x <- lapply(x, decode)
        y <- list(
            nCount = 46839L,
            nFeature = 245L,
            nMito = NA_integer_
        )
        expect_identical(x, y)
    }
})

test_that("SummarizedExperiment", {
    for (object in list(
        SummarizedExperiment = rse,
        SingleCellExperiment = sce
    )) {
        x <- calculateMetrics(object)
        expect_s4_class(x, "SummarizedExperiment")
    }
})

test_that("Low pass prefiltering", {
    ## All barcodes in example should pass.
    object <- sce
    x <- calculateMetrics(object, prefilter = TRUE)
    expect_identical(ncol(x), ncol(object))
    ## Simulate some poor barcodes.
    counts(object)[, seq_len(2L)] <- 0L
    x <- calculateMetrics(object, prefilter = TRUE)
    expect_identical(ncol(x), ncol(object) - 2L)
})
