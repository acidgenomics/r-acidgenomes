context("printString")

test_that("printString", {
    expect_identical(
        object = printString(c("hello", "world")),
        expected = "[1] \"hello\" \"world\""
    )
})
