## Multi-mapped identifiers that we can resolve:
## 2 16 19 22 23 25 28 29 31 34

## FIXME 9410 match failures: "11", "17", "62", "63", "64", "66",
"67", "68", "69", "73"....

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
