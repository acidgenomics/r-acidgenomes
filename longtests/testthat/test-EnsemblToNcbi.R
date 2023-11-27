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

## FIXME 151 mismatches between HGNC and Ensembl currently:
## - ENSG00000111215
## - ENSG00000169627
## - ENSG00000170667
## ...

test_that("Mapping consistency between Ensembl and HGNC", {
    hgnc <- Hgnc()
    ens <- makeGRangesFromEnsembl(organism = "Homo sapiens")
    x <- EnsemblToNcbi(hgnc)
    x <- as(x, "DFrame")
    metadata(x) <- list()
    rownames(x) <- NULL
    y <- EnsemblToNcbi(ens)
    y <- as(y, "DFrame")
    metadata(y) <- list()
    rownames(y) <- NULL
    genes <- sort(intersect(x[[1L]], y[[1L]]))
    x <- x[match(genes, table = x[[1L]]), ]
    y <- y[match(genes, table = y[[1L]]), ]
    expect_identical(x, y)
    ## > y[which(x[[2]] != y[[2]]), ]
})
