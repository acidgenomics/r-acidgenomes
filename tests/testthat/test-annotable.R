context("annotable")

skip_if_not(hasInternet())

test_that("genes", {
    object <- annotable(
        organism = organism,
        release = release,
        level = "genes"
    )
    expect_s3_class(object, "tbl_df")
    expect_identical(dim(object), c(67667L, 13L))
    expect_identical(object[["geneId"]][[1L]], "ENSG00000000003")
})
