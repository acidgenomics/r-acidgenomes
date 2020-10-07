context("selectSamples")

test_that("SummarizedExperiment", {
    object <- selectSamples(rse, condition = "A")
    expect_s4_class(object, class = "SummarizedExperiment")
    expect_identical(
        object = colnames(object),
        expected = paste0("sample0", seq_len(6L))
    )
})

test_that("SingleCellExperiment", {
    object <- selectSamples(sce, sampleID = "sample1")
    expect_identical(
        object = sampleNames(object),
        expected = c(sample1 = "sample1")
    )
    expect_identical(
        object = rownames(sampleData(object)),
        expected = "sample1"
    )
})
