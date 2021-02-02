context("organism")

## nolint start
GRanges <- gr
DataFrame <- as(as.data.frame(GRanges), "DataFrame")
matrix <- as.matrix(DataFrame)
## nolint end

test_that("organism", {
    for (object in list(
        "matrix" = matrix,
        "DataFrame" = DataFrame,
        "GRanges" = GRanges
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
