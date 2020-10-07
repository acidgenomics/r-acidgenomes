context("mapCellsToSamples")

test_that("2 samples", {
    samples <- paste0("sample", seq_len(2L))
    cells <- paste(samples, c("AAAAAAAA", "CCCCCCCC"), sep = "_")
    x <- mapCellsToSamples(cells = cells, samples = samples)
    expect_is(x, "factor")
    expect_identical(
        object = x,
        expected = as.factor(c(
            sample1_AAAAAAAA = "sample1",
            sample2_CCCCCCCC = "sample2"
        ))
    )
})

## Support for this is useful when working with example Seurat objects.
## Consider making this stricter and checking for ACGT specifically.
test_that("1 sample", {
    expect_identical(
        object = mapCellsToSamples(
            cells = c("AAAAAAAA", "CCCCCCCC"),
            samples = "sample"
        ),
        expected = as.factor(c(
            AAAAAAAA = "sample",
            CCCCCCCC = "sample"
        ))
    )
})

test_that("Invalid barcodes", {
    expect_error(
        object = mapCellsToSamples(cells = c("aaa", "bbb"), samples = "sample"),
        regexp = "ACGT"
    )
})

test_that("No separator", {
    samples <- paste0("sample", seq_len(2L))
    cells <- c("AAAAAAAA", "CCCCCCCC")
    expect_error(
        object = mapCellsToSamples(cells = cells, samples = samples),
        regexp = "allAreMatchingRegex"
    )
})

test_that("Match failure", {
    samples <- paste0("sample", seq_len(2L))
    cells <- paste("xxx", c("AAAAAAAA", "CCCCCCCC"), sep = "_")
    expect_error(
        object = mapCellsToSamples(cells = cells, samples = samples),
        expected = "\"sample1\" sample failed to match any cells."
    )
})
