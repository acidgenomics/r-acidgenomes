## This code depends on biomaRt and Ensembl, which can time out.

context("mapHumanOrthologs")

## FIXME This isn't sanitizing NA values the way we'd expect.

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
    if (is(object, "error")) {
        msg <- as.character(object)
        skip_if(
            condition = grepl(pattern = "biomaRt", x = msg),
            message = msg
        )
    }
    expected <- DataFrame(
        "geneId" = c(
            "ENSMUSG00000000001",
            "ENSMUSG00000000028",
            "ENSMUSG00000000037",
            "ENSMUSG00000000049"
        ),
        "geneName" = c(
            "Gnai3",
            "Cdc45",
            "Scml2",
            "Apoh"
        ),
        "humanGeneId" = c(
            "ENSG00000065135",
            "ENSG00000093009",
            "ENSG00000102098",
            "ENSG00000091583"
        ),
        "humanGeneName" = c(
            "GNAI3",
            "CDC45",
            "SCML2",
            "APOH"
        ),
        row.names = c(
            "ENSMUSG00000000001",
            "ENSMUSG00000000028",
            "ENSMUSG00000000037",
            "ENSMUSG00000000049"
        )
    )
    expect_identical(object, expected)
})
