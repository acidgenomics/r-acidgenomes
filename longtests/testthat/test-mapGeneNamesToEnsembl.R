hgnc <- Hgnc()
keep <- !is.na(hgnc[["ensemblGeneId"]])
hgnc <- hgnc[keep, ]
hgncGeneNames <- unlist(
    x = list(
        hgnc[["geneName"]],
        unlist(hgnc[["aliasSymbol"]]),
        unlist(hgnc[["prevSymbol"]])
    ),
    recursive = TRUE,
    use.names = FALSE
)
hgncGeneNames <- sort(unique(hgncGeneNames))

ncbi <- NcbiGeneInfo(organism = "Homo sapiens")
keep <- any(startsWith(x = ncbi[["dbXrefs"]], prefix = "Ensembl:"))
ncbi <- ncbi[keep, ]
dbXrefs <- ncbi[["dbXrefs"]]
dbXrefs <- dbXrefs[startsWith(x = dbXrefs, prefix = "Ensembl:")]
keep <- lengths(dbXrefs) == 1L
ncbi <- ncbi[keep, ]
ncbiGeneNames <- unlist(
    x = list(
        ncbi[["geneName"]],
        unlist(ncbi[["geneSynonyms"]])
    ),
    recursive = TRUE,
    use.names = FALSE
)
ncbiGeneNames <- sort(unique(ncbiGeneNames))

test_that("HGNC gene names", {
    expect_no_error(
        object = mapGeneNamesToEnsembl(
            genes = hgncGeneNames,
            organism = "Homo sapiens",
            hgnc = hgnc
        )
    )
    expect_error(
        object = mapGeneNamesToEnsembl(
            genes = hgncGeneNames,
            organism = "Homo sapiens",
            ncbi = ncbi
        ),
        regexp = "mapping failure"
    )
})

test_that("NCBI gene names", {
    expect_no_error(
        object = mapGeneNamesToEnsembl(
            genes = ncbiGeneNames,
            organism = "Homo sapiens",
            ncbi = ncbi
        )
    )
    expect_error(
        object = mapGeneNamesToEnsembl(
            genes = ncbiGeneNames,
            organism = "Homo sapiens",
            hgnc = hgnc
        ),
        regexp = "mapping failure"
    )
})

test_that("Case sensitivity", {
    expect_identical(
        object = mapGeneNamesToEnsembl(
            genes = c("HDGFL3", "Hdgfrp3"),
            organism = "Homo sapiens",
            ignoreCase = TRUE,
            hgnc = hgnc
        ),
        expected = rep("ENSG00000166503", times = 2L)
    )
    expect_error(
        object = mapGeneNamesToEnsembl(
            genes = c("HDGFL3", "HDGFRP3"),
            organism = "Homo sapiens",
            ignoreCase = FALSE,
            hgnc = hgnc
        ),
        regexp = "HDGFRP3"
    )
})
