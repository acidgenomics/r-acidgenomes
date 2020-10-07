context("convertTranscriptsToGenes")

skip_if_not(hasInternet())

tx2gene <-
    makeTx2GeneFromEnsembl(organism = "Homo sapiens", release = 87L)

test_that("character", {
    expect_identical(
        object = convertTranscriptsToGenes(
            object = c("ENST00000000233", "ENST00000000412"),
            tx2gene = tx2gene
        ),
        expected = factor(c(
            ENST00000000233 = "ENSG00000004059",
            ENST00000000412 = "ENSG00000003056"
        ))
    )
})

test_that("matrix", {
    object <- matrix(
        data = seq_len(8L),
        byrow = TRUE,
        nrow = 4L,
        ncol = 2L,
        dimnames = list(
            c(
                "ENST00000373020",
                "ENST00000494424",
                "ENST00000373031",
                "ENST00000485971"
            ),
            c("sample1", "sample2")
        )
    )
    expected <- matrix(
        data = c(4L, 6L, 12L, 14L),
        byrow = TRUE,
        nrow = 2L,
        ncol = 2L,
        dimnames = list(
            c("ENSG00000000003", "ENSG00000000005"),
            c("sample1", "sample2")
        )
    )
    expect_identical(
        object = convertTranscriptsToGenes(
            object = object,
            tx2gene = tx2gene,
            aggregate = TRUE
        ),
        expected = expected
    )
})

test_that("Invalid params", {
    expect_error(
        object = convertTranscriptsToGenes(
            object = "ENST00000000000",
            tx2gene = tx2gene
        ),
        regexp = "Failed to match transcripts: ENST00000000000"
    )
    expect_error(
        object = convertTranscriptsToGenes(c("ENST00000000233", NA)),
        regexp = "isCharacter"
    )
    expect_error(
        object = convertTranscriptsToGenes(c("ENST00000000233", "")),
        regexp = "isCharacter"
    )
})

test_that("SummarizedExperiment", {
    se <- convertTranscriptsToGenes(txse)
    expect_s4_class(se, "SummarizedExperiment")
    expect_identical(
        object = counts(se),
        expected = matrix(
            data = c(
                15L, 18L, 21L, 24L,
                51L, 54L, 57L, 60L
            ),
            nrow = 2L,
            ncol = 4L,
            byrow = TRUE,
            dimnames = list(
                c("ENSG00000000003", "ENSG00000000419"),
                c("sample1", "sample2", "sample3", "sample4")
            )
        )
    )
})
