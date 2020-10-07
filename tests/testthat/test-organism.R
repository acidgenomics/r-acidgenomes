context("organism")

object <- rse
rownames(object) <- as.character(rowData(rse)[["geneID"]])

## nolint start
matrix <- assay(object)
GRanges <- rowRanges(object)
DataFrame <- as(as.data.frame(GRanges), "DataFrame")
## nolint end

test_that("organism", {
    for (object in list(
        matrix = matrix,
        DataFrame = DataFrame,
        GRanges = GRanges,
        SummarizedExperiment = object
    )) {
        expect_identical(
            object = organism(object),
            expected = "Homo sapiens"
        )
    }
})

test_that("SE metadata stash", {
    org <- "xxx"
    metadata(object)[["organism"]] <- org
    expect_identical(organism(object), org)
})

rm(object)
