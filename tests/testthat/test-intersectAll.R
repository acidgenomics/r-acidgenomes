context("intersectAll")

test_that("list", {
    object <- list(
        a = c("a", "b", "c", "d", "e", "f"),
        b = c("b", "c", "d", "e", "f", "g"),
        c = c("c", "d", "e", "f", "g", "h")
    )
    expect_identical(
        object = intersectAll(object),
        expected = c("c", "d", "e", "f")
    )
})
