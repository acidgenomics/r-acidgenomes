test_that("Complete cases", {
    df <- DataFrame(
        "txId" = c(
            "transcript1",
            "transcript2",
            "transcript3",
            "transcript4",
            NA_character_
        ),
        "geneId" = c(
            "gene1",
            "gene1",
            "gene2",
            NA_character_,
            "gene3"
        ),
        row.names = LETTERS[seq_len(5L)]
    )
    object <- TxToGene(df)
    expect_s4_class(object, "TxToGene")
    expected <- DataFrame(
        "txId" = c(
            "transcript1",
            "transcript2",
            "transcript3"
        ),
        "geneId" = c(
            "gene1",
            "gene1",
            "gene2"
        ),
        row.names = LETTERS[seq_len(3L)]
    )
    metadata(expected) <- list(
        "date" = Sys.Date(),
        "dropped" = c("D" = 4L, "E" = 5L),
        "packageVersion" = .pkgVersion
    )
    expected <- new(Class = "TxToGene", expected)
    expect_identical(object, expected)
    expect_output(
        object = summary(object),
        regexp = "genes: 2"
    )
})

test_that("Ensure transcripts return sorted", {
    df <- DataFrame(
        "txId" = c(
            "transcript2",
            "transcript1",
            "transcript4",
            "transcript3"
        ),
        "geneId" = c(
            "gene1",
            "gene2",
            "gene1",
            "gene2"
        ),
        row.names = LETTERS[seq_len(4L)]
    )
    object <- TxToGene(df)
    expected <- DataFrame(
        "txId" = c(
            "transcript1",
            "transcript2",
            "transcript3",
            "transcript4"
        ),
        "geneId" = c(
            "gene2",
            "gene1",
            "gene2",
            "gene1"
        ),
        row.names = c("B", "A", "D", "C")
    )
    metadata(expected) <- list(
        "date" = Sys.Date(),
        "packageVersion" = .pkgVersion
    )
    expected <- new(Class = "TxToGene", expected)
    expect_identical(object, expected)
})
