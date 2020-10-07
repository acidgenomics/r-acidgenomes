context("humanize")

test_that("SummarizedExperiment", {
    samples <- letters[seq_along(colnames(rse))]
    names(samples) <- colnames(rse)
    sampleNames(rse) <- samples
    x <- humanize(rse)
    expect_identical(
        object = lapply(dimnames(x), head),
        expected = list(
            c("TSPAN6", "TNMD", "DPM1", "SCYL3", "C1orf112", "FGR"),
            c("a", "b", "c", "d", "e", "f")
        )
    )
})

## Check for unmodified return on SCE objects.
test_that("SingleCellExperiment", {
    x <- humanize(sce)
    expect_s4_class(x, "SingleCellExperiment")
    expect_identical(
        dimnames(x),
        list(
            as.character(mcols(rowRanges(x))[["geneName"]]),
            colnames(sce)
        )
    )
})
