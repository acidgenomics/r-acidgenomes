test_that("HGNC2Ensembl", {
    object <- HGNC2Ensembl()
    expect_s4_class(object, "HGNC2Ensembl")
    expect_identical(
        object = as.data.frame(object)["5", , drop = FALSE],
        expected = data.frame(
            "hgncId" = 5L,
            "ensemblId" = "ENSG00000121410",
            row.names = "5"
        )
    )
})
