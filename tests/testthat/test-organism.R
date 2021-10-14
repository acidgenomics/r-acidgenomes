context("organism")

## nolint start
GRanges <- gr
DataFrame <- as(as.data.frame(GRanges), "DFrame")
matrix <- as.matrix(DataFrame)
## nolint end

test_that("organism", {
    for (object in list(
        "DFrame" = DFrame,
        "GRanges" = GRanges,
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
