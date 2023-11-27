## Entrez identifiers that multimap:
## 2 16 19 22 23 25 28 29 31 34

test_that("character : all HGNC genes", {
    hgnc <- Hgnc()
    object <- sort(unique(na.omit(hgnc[["ncbiGeneId"]])))
    expect_error(
        object = NcbiToEnsembl(
            object = object,
            organism = "Homo sapiens"
        ),
        regexp = "match failures"
    )
})
