test_that("MGI2Ensembl", {
    object <- MGI2Ensembl()
    expect_s4_class(object, "MGI2Ensembl")
    expect_identical(
        object = as.data.frame(object)["87853", , drop = FALSE],
        expected = data.frame(
            "mgiId" = 87853L,
            "ensemblId" = "ENSMUSG00000027596",
            row.names = "87853"
        )
    )
})
