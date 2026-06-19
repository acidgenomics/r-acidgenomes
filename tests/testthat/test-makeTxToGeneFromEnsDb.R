skip_if_not_installed("EnsDb.Hsapiens.v75")

test_that("TxToGene from EnsDb", {
    object <- TxToGene(
        makeGRangesFromEnsDb(
            object = "EnsDb.Hsapiens.v75",
            level = "transcripts",
            ignoreVersion = FALSE
        )
    )
    expect_s4_class(object, "TxToGene")
    expect_identical(nrow(object), 215170L)
})
