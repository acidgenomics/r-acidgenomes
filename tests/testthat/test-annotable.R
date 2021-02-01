context("annotable")

skip_if_not(hasInternet())

test_that("genes", {
    object <- annotable(
        organism = organism,
        release = ensemblRelease,
        level = "genes",
        ignoreVersion = TRUE
    )
    expect_s3_class(object, "tbl_df")
    expect_identical(dim(object), c(68001L, 14L))
    expect_identical(object[["geneId"]][[1L]], "ENSG00000228572")
})
