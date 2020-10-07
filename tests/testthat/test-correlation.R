context("correlation")

data(
    correlation,
    package = "AcidTest",
    envir = environment()
)

list <- correlation
rm(correlation)

score <- 0.964

test_that("vector", {
    x <- list[["vector_x"]]
    y <- list[["vector_y"]]
    object <- correlation(x = x, y = y)
    object <- round(object, digits = 4L)
    expect_identical(object, score)
})

test_that("matrix", {
    x <- list[["matrix_x"]]
    y <- list[["matrix_y"]]
    expect_identical(
        object = cor(x),
        expected = correlation(x)
    )
    expect_identical(
        object = cor(x = c(x), y = c(y)),
        expected = correlation(x = x, y = y)
    )
})

test_that("SummarizedExperiment", {
    x <- list[["SummarizedExperiment_x"]]
    y <- list[["SummarizedExperiment_y"]]
    expect_identical(
        object = correlation(x = x, i = 1L),
        expected = correlation(assay(x, i = 1L))
    )
    expect_identical(
        object = round(correlation(x = x, i = 1L, j = 2L), digits = 4L),
        expected = score
    )
    expect_identical(
        object = round(correlation(x = x, y = y), digits = 4L),
        expected = score
    )
})
