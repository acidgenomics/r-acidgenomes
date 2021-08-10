## FIXME Need to improve code coverage of NA value handling.



context("Tx2Gene")

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
    object <- Tx2Gene(df)
    expect_s4_class(object, "Tx2Gene")
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
    metadata(expected) <- list("dropped" = c("D" = 4L, "E" = 5L))
    expected <- new(Class = "Tx2Gene", expected)
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
    object <- Tx2Gene(df)
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
    expected <- new(Class = "Tx2Gene", expected)
    expect_identical(object, expected)
})
