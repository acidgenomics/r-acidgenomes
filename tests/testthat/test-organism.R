context("organism")

## nolint start
GRanges <- gr  # FIXME
DataFrame <- as(as.data.frame(GRanges), "DataFrame")  # FIXME
matrix <- as.matrix(DataFrame)
## nolint end

test_that("organism", {
    for (object in list(
        "DataFrame" = DFrame,  # FIXME
        "GenomicRanges" = GRanges,  # FIXME
        "matrix" = matrix
    )) {
        expect_identical(
            object = organism(object),
            expected = "Homo sapiens"
        )
    }
})

test_that("Metadata stash", {
    object <- gr
    org <- "xxx"
    metadata(object)[["organism"]] <- org
    expect_identical(organism(object), org)
})
