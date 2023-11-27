## FIXME Add coverage for versioned ensembl identifiers.

## FIXME Move this to main tests when finished.

## FIXME Need to figure out ensembl genes that have tricky multi-mapping
## that doesn't match Hgnc.
##
## FIXME Need to figure out ncbi genes that have tricky multi-mapping that
## doesn't match Hgnc.

## Ensembl identifiers that multi-map:
## [1] "ENSG00000004866" "ENSG00000063587" "ENSG00000065615"
## [4] "ENSG00000076554" "ENSG00000088298" "ENSG00000090857"

## FIXME How to handle these match failures:
## 7461 match failures: "ENSG00000018607",
## "ENSG00000042304", "ENSG00000064489", "ENSG00000067601",
## "ENSG00000073905", "ENSG00000078319", "ENSG00000099251",
## "ENSG00000100058", "ENSG00000100068",
## "ENSG00000100181"....

test_that("character : all HGNC genes", {
    hgnc <- Hgnc()
    object <- sort(unique(na.omit(hgnc[["ensemblGeneId"]])))
    expect_error(
        object = EnsemblToNcbi(
            object = object,
            organism = "Homo sapiens"
        ),
        regexp = "match failures"
    )
})
