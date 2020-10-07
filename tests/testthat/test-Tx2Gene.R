context("Tx2Gene")

object <- Tx2Gene(txse)

test_that("Tx2Gene", {
    expect_s4_class(object, "Tx2Gene")
    expect_identical(dim(object), c(6L, 2L))
})

test_that("summary", {
    expect_output(
        object = summary(object),
        regexp = "organism: Homo sapiens"
    )
})
