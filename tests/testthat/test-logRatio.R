context("logRatio")

fc <- c(-8L, -4L, -2L, 1L, 2L, 4L, 8L)
lr <- seq(-3L, 3L, 1L)
vec1 <- seq(from = 1L, to = 5L, by = 1L)
vec2 <- vec1 ^ 2L
means <- c(vec1 = 2.605171, vec2 = 6.786916)

test_that("foldChangeToLogRatio", {
    expect_equal(foldChangeToLogRatio(fc), lr)
    expect_error(
        object = foldChangeToLogRatio(lr, base = 0L),
        regexp = "isPositive"
    )
})

test_that("logRatioToFoldChange", {
    expect_equal(logRatioToFoldChange(lr), fc)
    expect_error(
        object = logRatioToFoldChange(lr, base = 0L),
        regexp = "isPositive"
    )
})
