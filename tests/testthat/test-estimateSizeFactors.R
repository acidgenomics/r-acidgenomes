context("estimateSizeFactors")

test_that("SummarizedExperiment", {
    object <- rse
    expect_null(sizeFactors(object))
    object <- estimateSizeFactors(object)
    expect_type(sizeFactors(object), "double")
    expect_identical(
        unname(sizeFactors(object)),
        colData(object)[["sizeFactor"]]
    )
})

test_that("SingleCellExperiment", {
    object <- sce
    expect_null(sizeFactors(object))
    object <- estimateSizeFactors(object)
    expect_type(sizeFactors(object), "double")
    ## Check that values are stashed in int_colData, not colData.
    ## This was changed again in Bioconductor 3.11 release.
    ## > expect_identical(
    ## >     unname(sizeFactors(object)),
    ## >     unname(object@int_colData[["size_factor"]])
    ## > )
    ## > expect_null(colData(object)[["sizeFactor"]])
})
