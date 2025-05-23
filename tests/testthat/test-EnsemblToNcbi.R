test_that("character : Homo sapiens", {
    genes <- c(
        "ENSG00000063587.15",
        "ENSG00000004866.22",
        "ENSG00000000005.6",
        "ENSG00000000003.16"
    )
    object <- EnsemblToNcbi(genes, organism = NULL)
    expected <- DataFrame(
        "ensemblGeneId" = c(
            "ENSG00000063587",
            "ENSG00000004866",
            "ENSG00000000005",
            "ENSG00000000003"
        ),
        "ncbiGeneId" = c(
            10838L,
            7982L,
            64102L,
            7105L
        ),
        row.names = genes
    )
    expect_identical(
        object = as.data.frame(object),
        expected = as.data.frame(expected)
    )
    expect_error(
        EnsemblToNcbi(
            object = c(
                "ENSG00000000003",
                "ENSG00000000004"
            ),
            organism = "Homo sapiens"
        ),
        regexp = "ENSG00000000004"
    )
    expect_error(
        object = EnsemblToNcbi(
            object = "ENSG00000000000",
            organism = "Homo sapiens"
        ),
        regexp = "match failure"
    )
})

test_that("character : Mus musculus", {
    genes <- c(
        "ENSG00000063587.15",
        "ENSG00000004866.22",
        "ENSG00000000005.6",
        "ENSG00000000003.16"
    )
    object <- EnsemblToNcbi(genes, organism = NULL)
    expected <- DataFrame(
        "ensemblGeneId" = c(
            "ENSG00000063587",
            "ENSG00000004866",
            "ENSG00000000005",
            "ENSG00000000003"
        ),
        "ncbiGeneId" = c(
            10838L,
            7982L,
            64102L,
            7105L
        ),
        row.names = genes
    )
    expect_identical(
        object = as.data.frame(object),
        expected = as.data.frame(expected)
    )
    expect_error(
        object = EnsemblToNcbi(
            object = "ENSG00000000000",
            organism = "Homo sapiens"
        ),
        regexp = "match failure"
    )
})

test_that("character : Caenorhabditis elegans", {
    genes <- c("WBGene00000898", "WBGene00004804")
    object <- EnsemblToNcbi(
        object = genes,
        organism = "Caenorhabditis elegans"
    )
    expected <- DataFrame(
        "ensemblGeneId" = genes,
        "ncbiGeneId" = c(175410L, 177343L),
        row.names = genes
    )
    expect_identical(
        object = as.data.frame(object),
        expected = as.data.frame(expected)
    )
    expect_error(
        object = EnsemblToNcbi(
            object = "WBGene00000000",
            organism = "Caenorhabditis elegans"
        ),
        regexp = "ENSEMBL"
    )
})

test_that("EnsemblGenes", {
    hs <- makeGRangesFromEnsembl(
        organism = "Homo sapiens",
        genomeBuild = "GRCh38",
        release = 110L,
        ignoreVersion = FALSE
    )
    x <- EnsemblToNcbi(hs, useCurated = TRUE)
    y <- EnsemblToNcbi(hs, useCurated = FALSE)
    expect_true(metadata(x)[["useCurated"]])
    expect_null(metadata(y)[["useCurated"]])
    expect_identical(nrow(x), 26606L)
    expect_identical(nrow(y), 26608L)
    expect_identical(
        object = setdiff(y[[1L]], x[[1L]]),
        expected = c("ENSG00000290723.1", "ENSG00000291109.1")
    )
    y <- y[rownames(x), ]
    idx <- which(x[[2L]] != y[[2L]])
    ## NOTE This value changes with rolling NCBI updates.
    expect_length(idx, 108L)
    expect_identical(
        object = head(data.frame(
            "ensemblGeneId" = x[["ensemblGeneId"]][idx],
            "ncbiGeneId1" = x[["ncbiGeneId"]][idx],
            "ncbiGeneId2" = y[["ncbiGeneId"]][idx]
        )),
        expected = data.frame(
            "ensemblGeneId" = c(
                "ENSG00000111215.12",
                "ENSG00000169627.9",
                "ENSG00000176797.4",
                "ENSG00000177693.5",
                "ENSG00000178934.5",
                "ENSG00000180525.14"
            ),
            "ncbiGeneId1" = c(
                11272L,
                654483L,
                414325L,
                26682L,
                653499L,
                414235L
            ),
            "ncbiGeneId2" = c(
                5554L,
                552900L,
                55894L,
                124906935L,
                3963L,
                124900286L
            )
        )
    )
})

test_that("Hgnc", {
    hgnc <- Hgnc()
    object <- EnsemblToNcbi(hgnc)
    expect_identical(
        object = as.data.frame(object[1L:5L, ]),
        expected = data.frame(
            "ensemblGeneId" = c(
                "ENSG00000000003",
                "ENSG00000000005",
                "ENSG00000000419",
                "ENSG00000000457",
                "ENSG00000000460"
            ),
            "ncbiGeneId" = c(
                7105L,
                64102L,
                8813L,
                57147L,
                55732L
            ),
            row.names = as.character(c(
                11858L,
                17757L,
                3005L,
                19285L,
                25565L
            ))
        )
    )
})

test_that("Mgi", {
    mgi <- Mgi()
    object <- EnsemblToNcbi(mgi)
    expect_identical(
        object = as.data.frame(object[1L:5L, ]),
        expected = data.frame(
            "ensemblGeneId" = c(
                "ENSMUSG00000000001",
                "ENSMUSG00000000003",
                "ENSMUSG00000000028",
                "ENSMUSG00000000031",
                "ENSMUSG00000000037"
            ),
            "ncbiGeneId" = c(
                14679L,
                54192L,
                12544L,
                14955L,
                107815L
            ),
            row.names = as.character(c(
                95773L,
                1860484L,
                1338073L,
                95891L,
                1340042L
            ))
        )
    )
})
