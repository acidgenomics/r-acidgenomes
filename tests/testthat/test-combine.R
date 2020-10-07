context("combine")

x <- rse
colnames(x) <- paste0(
    "sample",
    str_pad(
        string = seq_len(ncol(x)),
        width = 2L,
        pad = "0"
    )
)

y <- x
colnames(y) <- paste0(
    "sample",
    str_pad(
        string = seq(from = ncol(y) + 1L, to = ncol(y) * 2L),
        width = 2L,
        pad = "0"
    )
)

## Note that this works on rowRanges internally.
test_that("RangedSummarizedExperiment", {
    c <- combine(x, y)
    expect_s4_class(c, "RangedSummarizedExperiment")

    samples <- paste0(
        "sample",
        str_pad(
            string = seq_len(24L),
            width = 2L,
            side = "left",
            pad = "0"
        )
    )
    expect_identical(colnames(c), samples)

    colData <- DataFrame(
        condition = as.factor(rep(rep(c("A", "B"), each = 6L), times = 2L)),
        row.names = samples
    )
    ## This will error out due to elementMetadata difference otherwise.
    expect_identical(
        object = as.data.frame(colData(c)),
        expected = as.data.frame(colData)
    )
})

## Note that this works on rowData internally.
test_that("SummarizedExperiment", {
    x <- as(x, "SummarizedExperiment")
    expect_true("geneName" %in% colnames(rowData(x)))

    y <- as(y, "SummarizedExperiment")
    expect_true("geneName" %in% colnames(rowData(y)))

    c <- combine(x, y)
    expect_s4_class(c, "SummarizedExperiment")
})

test_that("Column data mismatches", {
    colData(x)[["batch"]] <- as.factor(rep(c("a", "b"), times = ncol(x) / 2L))
    colData(y)[["group"]] <- as.factor(rep(c("c", "d"), each = ncol(x) / 2L))
    c <- combine(x, y)
    expect_s4_class(c, "RangedSummarizedExperiment")
    expect_true(anyNA(colData(c)[["batch"]]))
})

test_that("Disjoint metadata", {
    metadata(x)[["test"]] <- TRUE
    expect_true("test" %in% names(metadata(x)))
    c <- combine(x, y)
    expect_s4_class(c, "RangedSummarizedExperiment")
    expect_false("test" %in% names(metadata(c)))
})

test_that("Non-identical metadata", {
    metadata(x)[["test"]] <- TRUE
    metadata(y)[["test"]] <- FALSE
    expect_true("test" %in% names(metadata(x)))
    c <- combine(x, y)
    expect_s4_class(c, "RangedSummarizedExperiment")
    expect_false("test" %in% names(metadata(c)))
})

rm(x, y)

test_that("SingleCellExperiment", {
    x <- sce
    colnames(x) <- paste0(
        "cell",
        str_pad(
            string = seq_len(ncol(x)),
            width = 4L,
            pad = "0"
        )
    )

    y <- x
    colnames(y) <- paste0(
        "cell",
        str_pad(
            string = seq_len(ncol(y)) + ncol(y),
            width = 4L,
            pad = "0"
        )
    )

    ## Increase the sample ID numbers.
    sampleID <- y[["sampleID"]]
    sampleID <- gsub("1$", "3", sampleID)
    sampleID <- gsub("2$", "4", sampleID)
    y[["sampleID"]] <- as.factor(sampleID)

    c <- combine(x, y)
    expect_s4_class(c, "SingleCellExperiment")
    expect_identical(
        sort(rownames(sampleData(c))),
        paste0("sample", seq_len(4L))
    )
})
