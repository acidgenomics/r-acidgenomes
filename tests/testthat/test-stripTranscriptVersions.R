context("stripTranscriptVersions")

test_that("character", {
    ## Return unmodified if not Ensembl transcript (ENS*T).
    ## For example, check and make sure *C. elegans* transcripts are preserved.
    expect_identical(
        object = stripTranscriptVersions(c(
            "ENSMUST00000119854.7",
            "EGFP.1",
            "2L52.1a",
            "2L52.2"
        )),
        expected = c(
            "ENSMUST00000119854",
            "EGFP.1",
            "2L52.1a",
            "2L52.2"
        )
    )
})

test_that("matrix", {
    object <- mat
    rownames(object) <- c(
        "ENSMUST00000000001.1",
        "ENSMUST00000000001.2",
        "ENSMUST00000000002.1",
        "EGFP.1"
    )
    expect_identical(
        object = object %>%
            stripTranscriptVersions() %>%
            rownames(),
        expected = c(
            "ENSMUST00000000001",
            "ENSMUST00000000001",  # Dupes allowed in matrix.
            "ENSMUST00000000002",
            "EGFP.1"
        )
    )
})
