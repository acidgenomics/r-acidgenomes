context("cell2sample")

test_that("SCE factor return", {
    x <- cell2sample(sce, return = "factor")
    expect_is(x, "factor")
    expect_identical(levels(x), c("sample1", "sample2"))
    expect_true(hasNames(x))
})

test_that("SCE tibble return", {
    x <- cell2sample(sce, return = "tbl_df")
    expect_s3_class(x, "tbl_df")
    expect_identical(colnames(x), c("cellID", "sampleID"))
})
