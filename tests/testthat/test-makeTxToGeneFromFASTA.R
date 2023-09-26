test_that("Ensembl", {
    file <- txFastas[["ensembl"]]
    object <- makeTxToGeneFromFASTA(
        file = file,
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "TxToGene")
    expect_identical(nrow(object), 205131L)
    expect_identical(
        object = as.data.frame(object)[seq_len(3L), ],
        expected = data.frame(
            "txId" = c(
                "ENST00000000233.10",
                "ENST00000000412.8",
                "ENST00000000442.11"
            ),
            "geneId" = c(
                "ENSG00000004059.11",
                "ENSG00000003056.8",
                "ENSG00000173153.17"
            )
        )
    )
    object <- makeTxToGeneFromFASTA(
        file = file,
        ignoreVersion = TRUE
    )
    expect_identical(
        object = as.data.frame(object)[seq_len(3L), ],
        expected = data.frame(
            "txId" = c(
                "ENST00000000233",
                "ENST00000000412",
                "ENST00000000442"
            ),
            "geneId" = c(
                "ENSG00000004059",
                "ENSG00000003056",
                "ENSG00000173153"
            )
        )
    )
})

test_that("FlyBase", {
    file <- txFastas[["flybase"]]
    object <- makeTxToGeneFromFASTA(file)
    expect_s4_class(object, "TxToGene")
    expect_identical(nrow(object), 30799L)
    expect_identical(
        object = as.data.frame(object)[seq_len(3L), ],
        expected = data.frame(
            "txId" = c(
                "FBtr0005088",
                "FBtr0006151",
                "FBtr0070000"
            ),
            "geneId" = c(
                "FBgn0260439",
                "FBgn0000056",
                "FBgn0031081"
            )
        )
    )
})

test_that("GENCODE", {
    file <- txFastas[["gencode"]]
    object <- makeTxToGeneFromFASTA(
        file = file,
        ignoreVersion = FALSE
    )
    expect_s4_class(object, "TxToGene")
    expect_identical(nrow(object), 252416L)
    expect_identical(
        object = as.data.frame(object)[seq_len(3L), ],
        expected = data.frame(
            "txId" = c(
                "ENST00000000233.10",
                "ENST00000000412.8",
                "ENST00000000442.11"
            ),
            "geneId" = c(
                "ENSG00000004059.11",
                "ENSG00000003056.8",
                "ENSG00000173153.17"
            )
        )
    )
    object <- makeTxToGeneFromFASTA(
        file = file,
        ignoreVersion = TRUE
    )
    expect_identical(
        object = as.data.frame(object)[seq_len(3L), ],
        expected = data.frame(
            "txId" = c(
                "ENST00000000233",
                "ENST00000000412",
                "ENST00000000442"
            ),
            "geneId" = c(
                "ENSG00000004059",
                "ENSG00000003056",
                "ENSG00000173153"
            )
        )
    )
})

test_that("WormBase", {
    file <- txFastas[["wormbase"]]
    object <- makeTxToGeneFromFASTA(file)
    expect_s4_class(object, "TxToGene")
    expect_identical(nrow(object), 31998L)
    expect_identical(
        object = as.data.frame(object)[seq_len(3L), ],
        expected = data.frame(
            "txId" = c(
                "2L52.1a.1",
                "2L52.1b.1",
                "2RSSE.1a.1"
            ),
            "geneId" = c(
                "WBGene00007063",
                "WBGene00007063",
                "WBGene00007064"
            )
        )
    )
})
