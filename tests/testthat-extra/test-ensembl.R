context("Ensembl annotations")

mat <- .ucscMatrix
keep <- as.logical(mat[, "ensembldb"])
mat <- mat[keep, , drop = FALSE]

## Note that we're running at transcript level here to check the gene merge.
##
## Potential issues:
## - *Caenorhabditis elegans*
##     - Invalid transcript IDs.
## - *Canis familiaris*
##     - Using the full *Canis lupus familiaris* won't match.
## - *Saccharomyces cerevisiae*
##     - Invalid gene IDs.
##     - Iinvalid transcript IDs.
##     - No gene names.

test_that("UCSC genome build remaps", {
    mapply(
        organism = mat[, "organism"],
        genomeBuild = mat[, "ucsc"],
        FUN = function(organism, genomeBuild) {
            object <- makeGRangesFromEnsembl(
                organism = organism,
                genomeBuild = genomeBuild,
                level = "transcripts"
            )
            expect_s4_class(object, "GRanges")
        },
        SIMPLIFY = FALSE
    )
})
