skip_if_not_installed("EnsDb.Hsapiens.v75")

test_that("makeTxToGeneFromEnsDb", {
    object <- makeTxToGeneFromEnsDb("EnsDb.Hsapiens.v75")
    expect_s4_class(object, "TxToGene")
    # FIXME Check the number of rows here.
})
