context("makeSampleData")

test_that("makeSampleData", {
    object <- data.frame(
        genotype = rep(c("control", "wildtype"), times = 2L),
        treatment = rep(c("vector", "RNAi"), each = 2L),
        sampleName = paste("sample", seq_len(4L)),
        row.names = paste0("GSM000000", seq_len(4L))
    )
    x <- makeSampleData(object)
    expect_s4_class(x, "DataFrame")
    expect_true(all(bapply(x, is.factor)))
})
