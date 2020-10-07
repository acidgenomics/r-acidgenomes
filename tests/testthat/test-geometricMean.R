context("geometricMean")

int <- seq(from = 1L, to = 5L, by = 1L)
num <- int ^ 2L
df <- data.frame(int, num)
mat <- as.matrix(df)
mean <- c(int = 2.605171, num = 6.786916)

test_that("integer", {
    expect_identical(
        round(geometricMean(int), digits = 6L),
        mean[["int"]]
    )
})

test_that("numeric", {
    expect_identical(
        round(geometricMean(num), digits = 6L),
        mean[["num"]]
    )
})

test_that("matrix", {
    expect_identical(
        round(geometricMean(mat), digits = 6L),
        mean
    )
})

test_that("NaN on negative numbers", {
    expect_identical(
        geometricMean(seq(from = -5L, to = 5L, by = 1L)),
        NaN
    )
})

test_that("Zero propagation", {
    expect_identical(
        geometricMean(
            seq(from = 0L, to = 5L, by = 1L),
            zeroPropagate = TRUE
        ),
        0L
    )
    expect_identical(
        round(geometricMean(
            seq(from = 1L, to = 5L, by = 1L),
            zeroPropagate = TRUE
        ), digits = 6L),
        2.605171
    )
})
