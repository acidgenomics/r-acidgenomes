context("tpm")

test_that("SummarizedExperiment", {
    se <- SummarizedExperiment(
        assays = list(
            tpm = matrix(
                data = seq_len(4L),
                nrow = 2L,
                ncol = 2L,
                byrow = TRUE
            )
        )
    )
    x <- tpm(se)
    expect_is(x, "matrix")
})
