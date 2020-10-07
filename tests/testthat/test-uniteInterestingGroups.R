context("uniteInterestingGroups")

object <- as(datasets::mtcars, "DataFrame")

test_that("Check that 'interestingGroups' column gets defined as factor.", {
    x <- uniteInterestingGroups(
        object = object,
        interestingGroups = c("vs", "am", "gear")
    )
    expect_identical(
        levels(x[["interestingGroups"]]),
        c("0:0:3", "0:1:4", "0:1:5", "1:0:3", "1:0:4", "1:1:4", "1:1:5")
    )
})

test_that("Error on missing groups.", {
    expect_error(
        object = uniteInterestingGroups(
            object = object,
            interestingGroups = c("XXX", "YYY")
        ),
        regexp = "isSubset.*interestingGroups"
    )
})
