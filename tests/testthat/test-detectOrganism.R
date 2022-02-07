context("detectOrganism")

test_that("Homo sapiens", {
    for (object in c(
        "Homo sapiens",
        "hsapiens",
        "GRCh38",
        "hg38",
        "ENSG00000000001",
        "ENST00000000001"
    )) {
        expect_identical(
            object = detectOrganism(object),
            expected = "Homo sapiens"
        )
    }
})

test_that("Mus musculus", {
    for (object in c(
        "Mus musculus",
        "mmusculus",
        "GRCm38",
        "mm10",
        "ENSMUSG00000000001",
        "ENSMUST00000000001"
    )) {
        expect_identical(
            object = detectOrganism(object),
            expected = "Mus musculus"
        )
    }
})

test_that("Rattus norvegicus", {
    for (object in c(
        "Rattus norvegicus",
        "rnorvegicus",
        "Rnor_6.0",
        "rn5",
        "ENSRNOG00000000001",
        "ENSRNOT00000000001"
    )) {
        expect_identical(
            object = detectOrganism(object),
            expected = "Rattus norvegicus"
        )
    }
})

test_that("Danio rerio", {
    for (object in c(
        "Danio rerio",
        "drerio",
        "GRCz10",
        "danRer10",
        "ENSDARG00000000001",
        "ENSDART00000000001"
    )) {
        expect_identical(
            object = detectOrganism(object),
            expected = "Danio rerio"
        )
    }
})

test_that("Drosophila melanogaster", {
    for (object in c(
        "Drosophila melanogaster",
        "dmelanogaster",
        "BDGP6",
        "dm6",
        "FBgn0000001",
        "FBtr0000001"
    )) {
        expect_identical(
            object = detectOrganism(object),
            expected = "Drosophila melanogaster"
        )
    }
})

test_that("Caenorhabditis elegans", {
    for (object in c(
        "Caenorhabditis elegans",
        "celegans",
        "WBcel235",
        "ce11",
        "WBGene00000001"
    )) {
        expect_identical(
            object = detectOrganism(object),
            expected = "Caenorhabditis elegans"
        )
    }
})

test_that("Gallus gallus", {
    for (object in c(
        "Gallus gallus",
        "ggallus",
        "ENSGALG00000000001",
        "ENSGALT00000000001"
    )) {
        expect_identical(
            object = detectOrganism(object),
            expected = "Gallus gallus"
        )
    }
})

test_that("Ovis aries", {
    for (object in c(
        "Ovis aries",
        "oaries",
        "Oar_v3.1",
        "oviAri3",
        "ENSOARG00000000001",
        "ENSOART00000000001"
    )) {
        expect_identical(
            object = detectOrganism(object),
            expected = "Ovis aries"
        )
    }
})

## See also SGD website (https://www.yeastgenome.org/) for details.
test_that("Saccharomyces cerevisiae", {
    for (object in c(
        "Saccharomyces cerevisiae",
        "scerevisiae",
        "R64-1-1",
        "sacCer3",
        "YDL234C",
        "YDL236W",
        "YDL240C-A",
        "YDL247W-A"
    )) {
        expect_identical(
            object = detectOrganism(object),
            expected = "Saccharomyces cerevisiae"
        )
    }
})

test_that("Multiple organisms", {
    ## Function matches only the first genome.
    expect_identical(
        object = detectOrganism(c(
            "ENSG00000000001",
            "ENSG00000000002",
            "ENSMUSG00000000001",
            "ENSMUSG00000000002"
        )),
        expected = "Homo sapiens"
    )
})

test_that("Unsupported organism", {
    expect_error(
        object = detectOrganism("XXX"),
        regexp = "Failed to detect organism"
    )
})
