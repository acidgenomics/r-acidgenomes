context("toStringUnique")

test_that("toStringUnique", {
    expect_identical(
        toStringUnique(c("hello", "world", NA, "hello", "world", NA)),
        "hello, world"
    )
})
