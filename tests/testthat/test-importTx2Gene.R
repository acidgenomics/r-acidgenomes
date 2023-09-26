args <- list(
    "file" = file.path("cache", "txToGene.csv"),
    "organism" = "Homo sapiens",
    "genomeBuild" = "GRCh38",
    "release" = 100L
)

test_that("No version stripping", {
    object <- do.call(
        what = importTxToGene,
        args = append(
            x = args,
            values = list(
                "ignoreVersion" = c("tx" = FALSE, "gene" = FALSE)
            )
        )
    )
    expect_s4_class(object, "TxToGene")
    expect_identical(
        object = as.data.frame(object[1L, ]),
        expected = data.frame(
            "txId" = "ENST00000415118.1",
            "geneId" = "ENSG00000223997.1"
        )
    )
})

test_that("Strip transcript and gene versions", {
    object <- do.call(
        what = importTxToGene,
        args = append(
            x = args,
            values = list(
                "ignoreVersion" = c("tx" = TRUE, "gene" = TRUE)
            )
        )
    )
    expect_s4_class(object, "TxToGene")
    expect_identical(
        object = as.data.frame(object[1L, ]),
        expected = data.frame(
            "txId" = "ENST00000415118",
            "geneId" = "ENSG00000223997"
        )
    )
})
