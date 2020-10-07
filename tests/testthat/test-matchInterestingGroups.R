context("matchInterestingGroups")

test_that("matchInterestingGroups", {
    object <- rse
    expect_identical(
        object = matchInterestingGroups(object),
        expected = interestingGroups(object)
    )
    expect_identical(
        object = matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups(object)[[1L]]
        ),
        expected = interestingGroups(object)[[1L]]
    )
})
