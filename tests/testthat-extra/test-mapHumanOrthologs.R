## This code depends on biomaRt and Ensembl, which can time out.

context("extra | mapHumanOrthologs")

test_that("mapHumanOrthologs", {
    genes <- c(
        "ENSMUSG00000000001", "ENSMUSG00000000003",
        "ENSMUSG00000000028", "ENSMUSG00000000031",
        "ENSMUSG00000000037", "ENSMUSG00000000049"
    )
    ## This depends on biomaRt, and has a tendency to time out.
    object <- tryCatch(
        expr = mapHumanOrthologs(
            genes = genes,
            organism = NULL,
            ensemblRelease = 87L
        ),
        error = function(e) e
    )
    ## Skip if connection timed out.
    if (is(map, "error")) {
        msg <- as.character(map)
        skip_if(
            condition = grepl(pattern = "biomaRt", x = msg),
            message = msg
        )
    }
    expected <- DataFrame(
        "geneId" = c(
            "ENSMUSG00000000001",
            "ENSMUSG00000000003",
            "ENSMUSG00000000028",
            "ENSMUSG00000000031",
            "ENSMUSG00000000037",
            "ENSMUSG00000000049"
        ),
        "hgncId" = c(
            "ENSG00000065135",
            NA,
            "ENSG00000093009",
            NA,
            "ENSG00000102098",
            "ENSG00000091583"
        ),
        "geneName" = c(
            "Gnai3",
            "Pbsn",
            "Cdc45",
            "H19",
            "Scml2",
            "Apoh"
        ),
        "hgncName" = c(
            "GNAI3",
            NA,
            "CDC45",
            NA,
            "SCML2",
            "APOH"
        ),
        row.names = genes
    )
    expect_identical(object, expected)
})
