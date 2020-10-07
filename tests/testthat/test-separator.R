context("separator")

test_that("separator", {
    expect_identical(
        object = separator(times = 2L),
        expected = "=="
    )
})
