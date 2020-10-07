context("sanitizeSampleData")

test_that("sanitizeSampleData", {
    object <- DataFrame(
        genotype = rep(c("wt", "ko"), times = 2L),
        batch = rep(seq_len(2L), each = 2L),
        row.names = paste0("sample", seq_len(4L))
    )

    ## `sampleName` column is required.
    object[["sampleName"]] <- NULL
    expect_error(
        object = sanitizeSampleData(object),
        regexp = "sampleName"
    )

    ## And `sampleName` column can't contain duplicates.
    object[["sampleName"]] <- "XXX"
    expect_error(
        object = sanitizeSampleData(object),
        regexp = "hasNoDuplicates"
    )

    ## All columns should return factor.
    object[["sampleName"]] <- paste("sample", seq_len(nrow(object)))
    data <- sanitizeSampleData(object)
    expect_is(data, "DataFrame")
    expect_true(all(vapply(
        X = data,
        FUN = is.factor,
        FUN.VALUE = logical(1L)
    )))
})
