context("emptyRanges")

test_that("default", {
    object <- emptyRanges("XXX")
    expect_identical(
        object = levels(seqnames(object)),
        expected = "unknown"
    )
    expect_identical(
        object = ranges(object),
        expected = IRanges(
            start = 1L,
            end = 100L,
            names = "XXX"
        )
    )
})

test_that("mcolnames", {
    object <- emptyRanges(
        "EGFP",
        seqname = "transgene",
        mcolnames = c("geneId", "geneName")
    )
    expect_identical(
        object = names(mcols(object)),
        expected = c("geneId", "geneName")
    )
})
