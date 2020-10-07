context("sampleNames")

test_that("SummarizedExperiment", {
    object <- sampleNames(rse)
    expected <- as.character(sampleData(rse)[["sampleName"]])
    names(expected) <- colnames(rse)
    expect_identical(object, expected)
})

test_that("SE assignment", {
    samples <- letters[seq_len(ncol(rse))]
    names(samples) <- colnames(rse)
    sampleNames(rse) <- samples
    expect_s4_class(rse, "RangedSummarizedExperiment")
    expect_identical(sampleNames(rse), samples)
})

test_that("SE assignment failures", {
    expect_error(
        sampleNames(rse) <- "xxx",
        "The names of 'value' are NULL."
    )
    expect_error(
        sampleNames(rse) <- c(aaa = "bbb"),
        "areSetEqual"
    )
})

test_that("sampleName column isn't factor", {
    colData(rse)[["sampleName"]] <- "xxx"
    expect_error(
        object = sampleNames(rse),
        expected = paste(
            "sampleData() requires a 'sampleName'",
            "factor column in 'colData()'."
        )
    )
})

test_that("SingleCellExperiment", {
    expect_identical(
        object = sort(sampleNames(sce)),
        expected = c(
            sample1 = "sample1",
            sample2 = "sample2"
        )
    )
})

test_that("SCE assignment", {
    oldSamples <- sampleNames(sce)
    newSamples <- letters[seq_along(oldSamples)]
    names(newSamples) <- names(oldSamples)
    sampleNames(sce) <- newSamples
    expect_identical(
        object = sampleNames(sce),
        expected = c(sample1 = "a", sample2 = "b")
    )
})
