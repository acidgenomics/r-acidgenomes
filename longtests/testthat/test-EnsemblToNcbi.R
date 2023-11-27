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

## DataFrame with 111 rows and 2 columns
##       ensemblGeneId ncbiGeneId
##         <character>  <integer>
## 1   ENSG00000111215       5554
## 2   ENSG00000169627     552900
## 3   ENSG00000176797      55894
## 4   ENSG00000177693  124906935
## 5   ENSG00000178934       3963
## ...             ...        ...
## 107 ENSG00000280263     440888
## 108 ENSG00000280623  124900476
## 109 ENSG00000281880     440034
## 110 ENSG00000285219  100506207
## 111 ENSG00000285551  124902436

test_that("Mapping consistency between HGNC and Ensembl", {
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
