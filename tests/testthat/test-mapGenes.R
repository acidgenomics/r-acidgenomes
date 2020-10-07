g2s <- Gene2Symbol(rse)
geneIDs <- head(g2s[["geneID"]])
geneNames <- head(g2s[["geneName"]])
rownames <- head(rownames(rse))



context("mapGenesToRownames")

test_that("SummarizedExperiment", {
    expect_identical(
        object = mapGenesToRownames(rse, genes = rownames),
        expected = c(
            gene001 = "gene001",
            gene002 = "gene002",
            gene003 = "gene003",
            gene004 = "gene004",
            gene005 = "gene005",
            gene006 = "gene006"
        )
    )
    expect_identical(
        object = mapGenesToRownames(rse, genes = geneIDs),
        expected = c(
            ENSG00000000003 = "gene001",
            ENSG00000000005 = "gene002",
            ENSG00000000419 = "gene003",
            ENSG00000000457 = "gene004",
            ENSG00000000460 = "gene005",
            ENSG00000000938 = "gene006"
        )
    )
    expect_identical(
        object = mapGenesToRownames(rse, genes = geneNames),
        expected = c(
            TSPAN6 = "gene001",
            TNMD = "gene002",
            DPM1 = "gene003",
            SCYL3 = "gene004",
            C1orf112 = "gene005",
            FGR = "gene006"
        )
    )
})



context("mapGenesToIDs")

test_that("SummarizedExperiment", {
    expect_identical(
        object = mapGenesToIDs(rse, genes = rownames),
        expected = c(
            gene001 = "ENSG00000000003",
            gene002 = "ENSG00000000005",
            gene003 = "ENSG00000000419",
            gene004 = "ENSG00000000457",
            gene005 = "ENSG00000000460",
            gene006 = "ENSG00000000938"
        )
    )
    expect_identical(
        object = mapGenesToIDs(rse, genes = geneIDs),
        expected = c(
            ENSG00000000003 = "ENSG00000000003",
            ENSG00000000005 = "ENSG00000000005",
            ENSG00000000419 = "ENSG00000000419",
            ENSG00000000457 = "ENSG00000000457",
            ENSG00000000460 = "ENSG00000000460",
            ENSG00000000938 = "ENSG00000000938"
        )
    )
    expect_identical(
        object = mapGenesToIDs(rse, genes = geneNames),
        expected = c(
            TSPAN6 = "ENSG00000000003",
            TNMD = "ENSG00000000005",
            DPM1 = "ENSG00000000419",
            SCYL3 = "ENSG00000000457",
            C1orf112 = "ENSG00000000460",
            FGR = "ENSG00000000938"
        )
    )
})



context("mapGenesToSymbols")

test_that("SummarizedExperiment", {
    expect_identical(
        object = mapGenesToSymbols(rse, genes = rownames),
        expected = c(
            gene001 = "TSPAN6",
            gene002 = "TNMD",
            gene003 = "DPM1",
            gene004 = "SCYL3",
            gene005 = "C1orf112",
            gene006 = "FGR"
        )
    )
    expect_identical(
        object = mapGenesToSymbols(rse, genes = geneIDs),
        expected = c(
            ENSG00000000003 = "TSPAN6",
            ENSG00000000005 = "TNMD",
            ENSG00000000419 = "DPM1",
            ENSG00000000457 = "SCYL3",
            ENSG00000000460 = "C1orf112",
            ENSG00000000938 = "FGR"
        )
    )
})
