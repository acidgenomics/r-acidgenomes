context("organism")

test_that("organism", {
    list <- list()
    list[["GenomicRanges"]] <- gr
    list[["DataFrame"]] <- as.DataFrame(gr)
    list[["matrix"]] <- as.matrix(list[["DataFrame"]])
    for (object in list) {
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
