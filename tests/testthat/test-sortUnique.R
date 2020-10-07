context("sortUnique")

test_that("sortUnique", {
    expect_identical(
        sortUnique(c("milk", "eggs", "eggs", NA)),
        c("eggs", "milk", NA)
    )
})
