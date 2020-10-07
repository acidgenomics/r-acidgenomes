context("intersectionMatrix")

test_that("list", {
    object <- list(
        a = c("a", "b", "c", "d", "e", "f"),
        b = c("b", "c", "d", "e", "f", "g"),
        c = c("c", "d", "e", "f", "g", "h")
    )
    mat <- intersectionMatrix(object)
    expect_identical(
        object = dimnames(mat),
        expected = list(
            c("a", "b", "c", "d", "e", "f", "g", "h"),
            c("a", "b", "c")
        )
    )
    expect_identical(
        object = as.integer(rowSums(mat)),
        expected = c(1L, 2L, 3L, 3L, 3L, 3L, 2L, 1L)
    )
})
