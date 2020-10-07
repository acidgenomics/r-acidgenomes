context("interestingGroups")

test_that("SummarizedExperiment", {
    expect_identical(
        object = interestingGroups(rse),
        expected = "condition"
    )
})

test_that("Object with no metadata", {
    object <- rse
    metadata(object) <- list()
    expect_identical(interestingGroups(object), NULL)
})

test_that("Assignment method", {
    object <- rse
    intgroup <- interestingGroups(object)[[1L]]
    interestingGroups(object) <- intgroup
    expect_identical(
        object = interestingGroups(object),
        expected = intgroup
    )
    expect_error(
        object = interestingGroups(object) <- "XXX",
        regexp = "XXX"
    )
})
