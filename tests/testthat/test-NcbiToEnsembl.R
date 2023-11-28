test_that("character : Homo sapiens", {
    genes <- c(64102L, 10838L, 7982L, 7105L)
    object <- NcbiToEnsembl(genes, organism = "Homo sapiens")
    expect_identical(
        object = as.data.frame(object),
        expected = data.frame(
            "ncbiGeneId" = genes,
            "ensemblGeneId" = c(
                "ENSG00000000005",
                "ENSG00000063587",
                "ENSG00000004866",
                "ENSG00000000003"
            ),
            row.names = as.character(genes)
        )
    )
    expect_error(
        object = NcbiToEnsembl(
            object = 0L,
            organism = "Homo sapiens"
        ),
        regexp = "match failure"
    )
})

test_that("character : Mus musculus", {
    genes <- c(14679L, 54192L, 12544L, 14955L)
    object <- NcbiToEnsembl(genes, organism = "Mus musculus")
    expect_identical(
        object = as.data.frame(object),
        expected = data.frame(
            "ncbiGeneId" = genes,
            "ensemblGeneId" = c(
                "ENSMUSG00000000001",
                "ENSMUSG00000000003",
                "ENSMUSG00000000028",
                "ENSMUSG00000000031"
            ),
            row.names = as.character(genes)
        )
    )
    expect_error(
        object = NcbiToEnsembl(
            object = 0L,
            organism = "Mus musculus"
        ),
        regexp = "match failure"
    )
})

test_that("character : Caenorhabditis elegans", {
    genes <- c(175410L, 177343L)
    object <- NcbiToEnsembl(genes, organism = "Caenorhabditis elegans")
    expect_identical(
        object = as.data.frame(object),
        expected = data.frame(
            "ncbiGeneId" = genes,
            "ensemblGeneId" = c(
                "WBGene00000898",
                "WBGene00004804"
            ),
            row.names = as.character(genes)
        )
    )
    expect_error(
        object = NcbiToEnsembl(
            object = 0L,
            organism = "Caenorhabditis elegans"
        ),
        regexp = "ENTREZID"
    )
})

test_that("Hgnc", {
    hgnc <- Hgnc()
    object <- NcbiToEnsembl(hgnc)
    expect_identical(
        object = as.data.frame(object[1L:5L, ]),
        expected = data.frame(
            "ncbiGeneId" = c(1L, 2L, 3L, 9L, 10L),
            "ensemblGeneId" = c(
                "ENSG00000121410",
                "ENSG00000175899",
                "ENSG00000256069",
                "ENSG00000171428",
                "ENSG00000156006"
            ),
            row.names = as.character(c(5L, 7L, 8L, 7645L, 7646L))
        )
    )
})

test_that("Mgi", {
    mgi <- Mgi()
    object <- NcbiToEnsembl(mgi)
    expect_identical(
        object = as.data.frame(object[1L:5L, ]),
        expected = data.frame(
            "ncbiGeneId" = c(
                11287L,
                11298L,
                11302L,
                11303L,
                11304L
            ),
            "ensemblGeneId" = c(
                "ENSMUSG00000030359",
                "ENSMUSG00000020804",
                "ENSMUSG00000025375",
                "ENSMUSG00000015243",
                "ENSMUSG00000028125"
            ),
            row.names = as.character(c(
                87854L,
                1328365L,
                1197518L,
                99607L,
                109424L
            ))
        )
    )
})
