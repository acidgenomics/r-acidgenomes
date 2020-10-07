context("mapHumanOrthologs")

skip_if_not(hasInternet())

test_that("Mus musculus", {
    genes <- c(
        "ENSMUSG00000000001",
        "ENSMUSG00000000003",
        "ENSMUSG00000000028"
    )
    ## This depends on biomaRt, and has a tendency to time out.
    map <- tryCatch(
        expr = mapHumanOrthologs(genes, ensemblRelease = NULL),
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
    expect_s4_class(map, "DataFrame")
    expect_identical(
        object = map,
        expected = DataFrame(
            geneID = genes,
            geneName = c(
                "Gnai3",
                "Pbsn",
                "Cdc45"
            ),
            hgncID = c(
                "ENSG00000065135",
                NA,
                "ENSG00000093009"
            ),
            hgncName = c(
                "GNAI3",
                NA,
                "CDC45"
            ),
            row.names = genes
        )
    )
})
