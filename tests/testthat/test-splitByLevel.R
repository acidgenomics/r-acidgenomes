context("splitByLevel : DataFrame")

df <- DataFrame(
    compound = relevel(
        factor(
            c(
                rep(x = "dmso", each = 3L),
                rep(x = c("ts_0001", "ts_0002"), each = 9L)
            )
        ),
        ref = "dmso"
    ),
    concentration = factor(
        c(
            rep(x = 0L, each = 3L),
            rep(x = c(0.1, 1L, 10L), each = 3L, times = 2L)
        )
    ),
    replicate = factor(
        rep(
            seq(from = 1L, to = 3L),
            times = 7L
        )
    )
)

test_that("compound, no ref", {
    list <- splitByLevel(df, f = "compound", ref = FALSE)
    expect_identical(
        object = unique(as.character(list[[1L]][["compound"]])),
        expected = "dmso"
    )
    expect_identical(
        object = unique(as.character(list[[2L]][["compound"]])),
        expected = "ts_0001"
    )
    expect_identical(
        object = unique(as.character(list[[2L]][["concentration"]])),
        expected = c("0.1", "1", "10")
    )
})

test_that("compound, with ref", {
    list <- splitByLevel(df, f = "compound", ref = TRUE)
    expect_identical(
        object = unique(as.character(list[[1L]][["compound"]])),
        expected = c("dmso", "ts_0001")
    )
    expect_identical(
        object = unique(as.character(list[[2L]][["compound"]])),
        expected = c("dmso", "ts_0002")
    )
    expect_identical(
        object = unique(as.character(list[[1L]][["concentration"]])),
        expected = c("0", "0.1", "1", "10")
    )
})

test_that("concentration, no ref", {
    list <- splitByLevel(df, f = "concentration", ref = FALSE)
    expect_identical(
        object = unique(as.character(list[[1L]][["compound"]])),
        expected = "dmso"
    )
    expect_identical(
        object = unique(as.character(list[[2L]][["compound"]])),
        expected = c("ts_0001", "ts_0002")
    )
    expect_identical(
        object = unique(as.character(list[[2L]][["concentration"]])),
        expected = "0.1"
    )
})

test_that("concentration, with ref", {
    list <- splitByLevel(df, f = "concentration", ref = TRUE)
    expect_identical(
        object = unique(as.character(list[[1L]][["compound"]])),
        expected = c("dmso", "ts_0001", "ts_0002")
    )
    expect_identical(
        object = unique(as.character(list[[1L]][["concentration"]])),
        expected = c("0", "0.1")
    )
    expect_identical(
        object = unique(as.character(list[[2L]][["concentration"]])),
        expected = c("0", "1")
    )
})
