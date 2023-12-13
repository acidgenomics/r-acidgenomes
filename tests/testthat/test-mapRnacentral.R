test_that("Ensembl", {
    object <- mapRnacentral(
        organism = "Homo sapiens",
        database = "Ensembl"
    )
    expect_identical(
        object = as.data.frame(object[1L, ]),
        expected = data.frame(
            "rnacentralId" = "URS0000000055",
            "database" = "ENSEMBL",
            "txId" = "ENST00000585414",
            "taxId" = 9606L,
            "txBiotype" = "lncRNA",
            "geneId" = "ENSG00000226803.9"
        )
    )
})

test_that("GENCODE", {
    object <- mapRnacentral(
        organism = "Homo sapiens",
        database = "GENCODE"
    )
    expect_identical(
        object = as.data.frame(object[1L, ]),
        expected = data.frame(
            "rnacentralId" = "URS0000000055",
            "database" = "ENSEMBL_GENCODE",
            "txId" = "ENST00000585414",
            "taxId" = 9606L,
            "txBiotype" = "lncRNA",
            "geneId" = "ENSG00000226803.9"
        )
    )
})

test_that("HGNC", {
    object <- mapRnacentral(
        organism = "Homo sapiens",
        database = "HGNC"
    )
    expect_identical(
        object = as.data.frame(object[1L, ]),
        expected = data.frame(
            "rnacentralId" = "URS0000000A8C",
            "database" = "HGNC",
            "txId" = "HGNC:12655",
            "taxId" = 9606L,
            "txBiotype" = "vault_RNA",
            "geneId" = "VTRNA1-2"
        )
    )
})

test_that("MGI", {
    object <- mapRnacentral(
        organism = "Mus musculus",
        database = "MGI"
    )
    expect_identical(
        object = as.data.frame(object[1L, ]),
        expected = data.frame(
            "rnacentralId" = "URS0000007529",
            "database" = "MGI",
            "txId" = "MGI:102485",
            "taxId" = 10090L,
            "txBiotype" = "tRNA",
            "geneId" = "mt-Th"
        )
    )
})

test_that("RefSeq", {
    object <- mapRnacentral(
        organism = "Homo sapiens",
        database = "RefSeq"
    )
    expect_identical(
        object = as.data.frame(object[1L, ]),
        expected = data.frame(
            "rnacentralId" = "URS0000000A8C",
            "database" = "REFSEQ",
            "txId" = "NR_026704",
            "taxId" = 9606L,
            "txBiotype" = "vault_RNA",
            "geneId" = "VTRNA1-2"
        )
    )
})
