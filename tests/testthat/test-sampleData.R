context("sampleData : SummarizedExperiment")

test_that("Return", {
    object <- rse
    ## Check that `sampleName` and `interestingGroups` auto-populate.
    expect_identical(
        object = setdiff(
            x = colnames(sampleData(object)),
            y = colnames(colData(object))
        ),
        expected = c("sampleName", "interestingGroups")
    )
    x <- sampleData(object)[, colnames(colData(object)), drop = FALSE]
    y <- colData(object)
    expect_identical(
        object = as.data.frame(x),
        expected = as.data.frame(y)
    )
    ## Empty `colData` is supported.
    object <- rse
    colData(object) <- DataFrame(row.names = colnames(object))
    interestingGroups(object) <- NULL
    expect_silent(sampleData(object))
})

test_that("Empty return", {
    se <- SummarizedExperiment()
    expect_identical(sampleData(se), DataFrame())
})

test_that("Assignment", {
    object <- rse
    sampleData(object)[["test"]] <- as.factor(seq_len(ncol(object)))
    expect_is(sampleData(object)[["test"]], "factor")
})



context("sampleData : SingleCellExperiment")

## Note that this doesn't sort the row names automatically. Here we're doing
## this internally in the check, to make the example "sce" object more resistant
## to code-breaking updates.
test_that("Return", {
    samples <- paste0("sample", seq_len(2L))
    sd <- sampleData(sce)
    expect_identical(sort(rownames(sd)), samples)
    expect_identical(
        sd[samples, ],
        DataFrame(
            sampleName = as.factor(samples),
            interestingGroups = as.factor(samples),
            row.names = samples
        )
    )
})

test_that("Empty return", {
    sce <- SingleCellExperiment()
    expect_identical(sampleData(sce), DataFrame())
})

test_that("No sample info", {
    colData(sce) <- DataFrame(row.names = colnames(sce))
    expect_identical(
        sampleData(sce),
        DataFrame(
            sampleName = factor("unknown"),
            interestingGroups = factor("unknown"),
            row.names = "unknown"
        )
    )
})

test_that("Assignment", {
    sd <- sampleData(sce) %>% .[sort(rownames(.)), , drop = FALSE]
    batch <- as.factor(seq_len(nrow(sd)))
    sd[["batch"]] <- batch
    sampleData(sce) <- sd
    sd <- sampleData(sce) %>% .[sort(rownames(.)), , drop = FALSE]
    samples <- paste0("sample", seq(2L))
    expect_identical(
        object = sd,
        expected = DataFrame(
            sampleName = as.factor(samples),
            batch = batch,
            interestingGroups = as.factor(samples),
            row.names = samples
        )
    )
})
