context("importTx2Gene")

args <- list(
    file = file.path("cache", "tx2gene.csv"),
    organism = "Homo sapiens",
    genomeBuild = "GRCh38",
    ensemblRelease = 100L
)

test_that("No version stripping", {
    object <- do.call(
        what = importTx2Gene,
        args = c(
            args,
            ignoreVersion = FALSE,
            ignoreGeneVersion = FALSE
        )
    )
    expect_is(object, "Tx2Gene")
    expect_identical(
        object = as.data.frame(object[1L, ]),
        expected = data.frame(
            transcriptID = "ENST00000631435.1",
            geneID = "ENSG00000282253.1"
        )
    )
})

test_that("Strip transcript and gene versions", {
    object <- do.call(
        what = importTx2Gene,
        args = c(
            args,
            ignoreVersion = TRUE,
            ignoreGeneVersion = TRUE
        )
    )
    expect_is(object, "Tx2Gene")
    expect_identical(
        object = as.data.frame(object[1L, ]),
        expected = data.frame(
            transcriptID = "ENST00000631435",
            geneID = "ENSG00000282253"
        )
    )
})
